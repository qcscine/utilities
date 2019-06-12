/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_IRCOPTIMIZER_H_
#define UTILS_IRCOPTIMIZER_H_

#include "Utils/CalculatorBasics/PropertyList.h"
#include "Utils/CalculatorBasics/Results.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/GeometryOptimization/GeometryOptimizer.h"
#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include <Core/Interfaces/Calculator.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @brief The base class for all IRC optimizers.
 *
 * The main purpose of this base class is to hide the template parameter(s)
 * of the derived class(es).
 */
class IRCOptimizerBase {
 public:
  static constexpr const char* ircInitialStepSizeKey = "irc_initial_step_size";
  static constexpr const char* ircTransfromCoordinatesKey = "irc_transfrom_coordinates";

  /// @brief Default constructor.
  IRCOptimizerBase() = default;
  /// @brief Virtual default destructor.
  virtual ~IRCOptimizerBase() = default;
  /**
   * @brief The main functionality of the IRC optimizer.
   *
   * This function wraps the optimize functions of the underlying optimizer.
   *
   * @param atoms The AtomCollection (Geometry) to be optimized.
   * @return int  The final number of optimization cycles carried out.
   */
  virtual int optimize(AtomCollection& atoms, const Eigen::VectorXd& mode, bool forward = true) = 0;
  /**
   * @brief Function to apply the given settings to underlying classes.
   * @param settings The new settings.
   */
  virtual void setSettings(const Settings& settings) = 0;
  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  virtual Settings getSettings() const = 0;
  /**
   * @brief Add an observer function that will be triggered in each iteration.
   *
   * @param function A function to be executed in every loop of the optimization.
   *                 The function will have access to the current cycle count,
   *                 the current value and to a const reference of the current
   *                 parameters.
   */
  virtual void addObserver(std::function<void(const int&, const double&, const Eigen::VectorXd&)> function) = 0;
  /**
   * @brief Clear all existing observer functions.
   *
   * For optimization problems with very fast evaluations of the underlying function
   * the removal of all observers can increase performance as the observers are given
   * as std::functions and can not be added via templates.
   */
  virtual void clearObservers() = 0;
  /// @brief The size of the initial step along the chosen mode.
  double initialStepSize = 0.01;
  /**
   * @brief Switch to transform the coordinates from Cartesian into an internal space.
   *
   * The optimization will be carried out in the internal coordinate space possibly
   * accellerating convergence.
   */
  bool transformCoordinates = true;
};

/**
 * @brief Settings for an IRCOptimizer.
 *
 * Uses template arguments in order to automatically include the
 * settings of underlying objects into the given settings.
 *
 * @tparam OptimizerType The underlying Optimizer class.
 * @tparam ConvergenceCheckType The underlying ConvergenceCheck class.
 */
template<class OptimizerType, class ConvergenceCheckType>
class IRCOptimizerSettings : public Settings {
 public:
  /**
   * @brief Construct a new IRCOptimizerSettings object.
   *
   * Sets the default values of the settings to the current values set in the objects
   * given to the constructor.
   *
   * @param ircBase The IRC optimizer.
   * @param optimizer The optimizer.
   * @param check The convergence check criteria.
   */
  IRCOptimizerSettings(const IRCOptimizerBase& ircBase, const OptimizerType& optimizer, const ConvergenceCheckType& check)
    : Settings("IRCOptimizerSettings") {
    optimizer.addSettingsDescriptors(this->_fields);
    check.addSettingsDescriptors(this->_fields);

    UniversalSettings::BoolDescriptor irc_initial_step_size("The size of the initial step along the chosen mode.");
    irc_initial_step_size.setDefaultValue(ircBase.initialStepSize);
    this->_fields.push_back(IRCOptimizerBase::ircInitialStepSizeKey, irc_initial_step_size);

    UniversalSettings::BoolDescriptor irc_transfrom_coordinates(
        "Switch to transform the coordinates from Cartesian into an internal space.");
    irc_transfrom_coordinates.setDefaultValue(ircBase.transformCoordinates);
    this->_fields.push_back(IRCOptimizerBase::ircTransfromCoordinatesKey, irc_transfrom_coordinates);
    this->resetToDefaults();
  }
};

/**
 * @brief A version of the GeometryOptimizer that optimizes along an internal reaction coordinate (IRC).
 *
 * This optimizer mass-weights the actual gradient in order to optimize in the mass-weighted
 * coordinate system.
 *
 * @tparam OptimizerType Expects any of the Optimizer classes. Note that some special optimizers
 *                       may not yet be supported or may need additional specialization.
 */
template<class OptimizerType>
class IRCOptimizer : public IRCOptimizerBase {
 public:
  /**
   * @brief Construct a new IRCOptimizer object.
   * @param calculator The calculator to be used for the single point/gradient calculations.
   */
  IRCOptimizer(Core::Calculator& calculator) : _calculator(calculator){};

  /**
   * @brief See IRCOptimizerBase::optimize().
   *
   * @param atoms The AtomCollection (Geometry) to be optimized.
   * @param mode  The mode to follow in the IRC.
   * @param forward A boolean signaling to follow the mode forwards (true, current positions + mode)
   *                or backwards (false, current positions - mode)
   * @return int  The final number of optimization cycles carried out.
   */
  virtual int optimize(AtomCollection& atoms, const Eigen::VectorXd& mode, bool forward = true) final {
    // Configure Calculator
    _calculator.setStructure(atoms);
    _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
    // Collect masses
    Eigen::VectorXd masses = Eigen::VectorXd::Zero(atoms.size());
    const auto& elements = atoms.getElements();
    for (unsigned int i = 0; i < atoms.size(); i++) {
      masses[i] = ElementInfo::mass(elements[i]);
    }
    masses.array() /= masses.sum();
    // Transformation into internal basis
    Eigen::MatrixXd transformation;
    if (this->transformCoordinates) {
      transformation = Geometry::calculateRotTransFreeTransformMatrix(atoms.getPositions(), atoms.getElements());
    }
    // Define update function
    const unsigned int nAtoms = atoms.size();
    auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
      Utils::PositionCollection coordinates;
      if (this->transformCoordinates) {
        auto tmp = (transformation * parameters).eval();
        coordinates = Eigen::Map<const Utils::PositionCollection>(tmp.data(), nAtoms, 3);
      }
      else {
        coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
      }
      _calculator.modifyPositions(coordinates);
      Utils::Results results = _calculator.calculate("Geometry Optimization Cycle");
      value = results.getEnergy();
      auto gradientMatrix = results.getGradients();
      gradientMatrix.col(0).array() *= masses.array();
      gradientMatrix.col(1).array() *= masses.array();
      gradientMatrix.col(2).array() *= masses.array();
      if (this->transformCoordinates) {
        auto tmp = Eigen::Map<const Eigen::VectorXd>(gradientMatrix.data(), nAtoms * 3);
        gradients = (transformation.transpose() * tmp).eval();
      }
      else {
        gradients = Eigen::Map<const Eigen::VectorXd>(gradientMatrix.data(), nAtoms * 3);
      }
    };
    // Move initial positions along mode
    Eigen::VectorXd positions;
    if (this->transformCoordinates) {
      Eigen::VectorXd tmp = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), nAtoms * 3);
      tmp += (forward ? this->initialStepSize : -1.0 * this->initialStepSize) * (mode / mode.norm());
      positions = (transformation.transpose() * tmp).eval();
    }
    else {
      positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), nAtoms * 3);
      positions += (forward ? this->initialStepSize : -1.0 * this->initialStepSize) * (mode / mode.norm());
    }
    // Optimize
    auto cycles = optimizer.optimize(positions, update, check);
    // Update Atom collection and return
    Utils::PositionCollection coordinates;
    if (this->transformCoordinates) {
      auto tmp = (transformation * positions).eval();
      coordinates = Eigen::Map<const Utils::PositionCollection>(tmp.data(), nAtoms, 3);
    }
    else {
      coordinates = Eigen::Map<const Utils::PositionCollection>(positions.data(), nAtoms, 3);
    }
    atoms.setPositions(coordinates);
    return cycles;
  }

  /**
   * @brief Function to apply the given settings to underlying classes.
   * @param settings The new settings.
   */
  virtual void setSettings(const Settings& settings) override {
    check.applySettings(settings);
    optimizer.applySettings(settings);
    this->initialStepSize = settings.getBool(IRCOptimizerBase::ircInitialStepSizeKey);
    this->transformCoordinates = settings.getBool(IRCOptimizerBase::ircTransfromCoordinatesKey);
  };
  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  virtual Settings getSettings() const override {
    return IRCOptimizerSettings<OptimizerType, GradientBasedCheck>(*this, optimizer, check);
  };
  /**
   * @brief Add an observer function that will be triggered in each iteration.
   *
   * @param function A function to be executed in every loop of the optimization.
   *                 The function will have access to the current cycle count
   *                 the current value and to a const reference of the current
   *                 parameters.
   */
  virtual void addObserver(std::function<void(const int&, const double&, const Eigen::VectorXd&)> function) final {
    optimizer.addObserver(function);
  }
  /**
   * @brief Clear all existing observer functions.
   *
   * For optimization problems with very fast evaluations of the underlying function
   * the removal of all observers can increase performance as the observers are given
   * as std::functions and can not be added via templates.
   */
  virtual void clearObservers() final {
    optimizer.clearObservers();
  }
  /// @brief The underlying optimizer, public in order to change it's settings.
  OptimizerType optimizer;
  /// @brief The underlying convergence check, public in order to change it's settings.
  GradientBasedCheck check;

 private:
  Core::Calculator& _calculator;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_GEOMETRYOPTIMIZER_H_
