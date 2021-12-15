/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_IRCOPTIMIZER_H_
#define UTILS_IRCOPTIMIZER_H_

#include "Utils/CalculatorBasics/PropertyList.h"
#include "Utils/CalculatorBasics/Results.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Geometry/InternalCoordinates.h"
#include "Utils/GeometryOptimization/GeometryOptimizer.h"
#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/GradientBased/SteepestDescent.h"
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
class IrcOptimizerBase {
 public:
  static constexpr const char* ircInitialStepSizeKey = "irc_initial_step_size";
  static constexpr const char* ircCoordinateSystemKey = "irc_coordinate_system";

  /// @brief Default constructor.
  IrcOptimizerBase() = default;
  /// @brief Virtual default destructor.
  virtual ~IrcOptimizerBase() = default;
  /**
   * @brief The main functionality of the IRC optimizer.
   *
   * This function wraps the optimize functions of the underlying optimizer.
   *
   * @param atoms The AtomCollection (Geometry) to be optimized.
   * @param log The logger to which eventual output is written.
   * @param mode  The mode to follow in the IRC.
   * @param forward A boolean signaling to follow the mode forwards (true, current positions + mode)
   *                or backwards (false, current positions - mode)
   * @return int  The final number of optimization cycles carried out.
   */
  virtual int optimize(AtomCollection& atoms, Core::Log& log, const Eigen::VectorXd& mode, bool forward = true) = 0;
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
   * @brief Get the settings of the calculator used for the energy calculations during the optimization.
   * @return std::shared_ptr<Settings> The settings of the calculator.
   */
  virtual const std::shared_ptr<Settings> getCalculatorSettings() const = 0;
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
  /**
   * @brief The size of the initial step along the chosen mode.
   *
   * Maximum step along one coordinate of one atom.
   */
  double initialStepSize = 0.3;
  /**
   * @brief Set the coordinate system in which the optimization shall be performed
   *
   * The optimization can be carried out in the internal coordinate space or with removed translations and rotations
   * possibly accelerating convergence.
   */
  CoordinateSystem coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
  /**
   * @brief The underlying convergence check
   *
   * @return GradientBasedCheck the class holding all convergence thresholds.
   */
  virtual const GradientBasedCheck getConvergenceCheck() const = 0;
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
  IRCOptimizerSettings(const IrcOptimizerBase& ircBase, const OptimizerType& optimizer, const ConvergenceCheckType& check)
    : Settings("IRCOptimizerSettings") {
    optimizer.addSettingsDescriptors(this->_fields);
    check.addSettingsDescriptors(this->_fields);

    UniversalSettings::DoubleDescriptor irc_initial_step_size("The size of the initial step along the chosen mode.");
    irc_initial_step_size.setDefaultValue(ircBase.initialStepSize);
    this->_fields.push_back(IrcOptimizerBase::ircInitialStepSizeKey, irc_initial_step_size);

    UniversalSettings::OptionListDescriptor irc_coordinate_system("Set the coordinate system.");
    irc_coordinate_system.addOption("internal");
    irc_coordinate_system.addOption("cartesianWithoutRotTrans");
    irc_coordinate_system.addOption("cartesian");
    irc_coordinate_system.setDefaultOption(CoordinateSystemInterpreter::getStringFromCoordinateSystem(ircBase.coordinateSystem));
    this->_fields.push_back(IrcOptimizerBase::ircCoordinateSystemKey, irc_coordinate_system);

    this->resetToDefaults();
  }
};

/**
 * @brief A version of the GeometryOptimizer that optimizes along an internal reaction coordinate (IRC).
 *
 * This optimizer mass-weights the actual gradient. The optimization is still carried out in the cartesian
 * coordinate system using the mass-weighted gradients. The mode for the initial displacement is obtained
 * from a mass-weighted Hessian, however back-scaled to cartesian coordinates.
 *
 * @tparam OptimizerType Expects any of the Optimizer classes. Note that some special optimizers
 *                       may not yet be supported or may need additional specialization.
 */
template<class OptimizerType>
class IrcOptimizer : public IrcOptimizerBase {
 public:
  /**
   * @brief Construct a new IRCOptimizer object.
   * @param calculator The calculator to be used for the single point/gradient calculations.
   */
  IrcOptimizer(Core::Calculator& calculator) : _calculator(calculator){};

  /**
   * @brief See IRCOptimizerBase::optimize().
   *
   * @param atoms The AtomCollection (Geometry) to be optimized.
   * @param log The logger to which eventual output is written.
   * @param mode  The mode to follow in the IRC.
   * @param forward A boolean signaling to follow the mode forwards (true, current positions + mode)
   *                or backwards (false, current positions - mode)
   * @return int  The final number of optimization cycles carried out.
   */
  int optimize(AtomCollection& atoms, Core::Log& log, const Eigen::VectorXd& mode, bool forward = true) final {
    // Collect square root of masses
    Eigen::VectorXd masses = Eigen::VectorXd::Zero(atoms.size());
    const unsigned int nAtoms = atoms.size();
    const auto& elements = atoms.getElements();
    for (unsigned int i = 0; i < nAtoms; i++) {
      masses[i] = sqrt(ElementInfo::mass(elements[i]));
    }
    // Distort the structure by the initial step
    // Get mode and scale it
    double maxVal = mode.array().abs().maxCoeff();
    Eigen::VectorXd modev = (forward ? (this->initialStepSize) : (-1.0 * this->initialStepSize)) * (mode / maxVal);
    auto modem = Eigen::Map<Utils::PositionCollection>(modev.data(), nAtoms, 3);
    // Move initial positions along mode
    Utils::PositionCollection coordinates = atoms.getPositions() + modem;
    atoms.setPositions(coordinates);
    // Configure Calculator
    _calculator.setStructure(atoms);
    _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
    // Transform into internal coordinates
    std::shared_ptr<InternalCoordinates> transformation = nullptr;
    if (this->coordinateSystem == CoordinateSystem::Internal) {
      transformation = std::make_shared<InternalCoordinates>(atoms);
    }
    else if (this->coordinateSystem == CoordinateSystem::CartesianWithoutRotTrans) {
      transformation = std::make_shared<InternalCoordinates>(atoms, true);
    }
    // Count cycles for eventual restart
    int cycle = 0;
    // Define update function
    auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
      cycle++;
      Utils::PositionCollection coordinates;
      if (transformation) {
        coordinates = transformation->coordinatesToCartesian(parameters);
      }
      else {
        coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
      }
      _calculator.modifyPositions(coordinates);
      Results results =
          CalculationRoutines::calculateWithCatch(_calculator, log, "Aborting optimization due to failed calculation");
      value = results.get<Property::Energy>();
      auto gradientMatrix = results.get<Property::Gradients>();
      // Mass weight the gradients by dividing with the square root of masses
      gradientMatrix.col(0).array() /= masses.array();
      gradientMatrix.col(1).array() /= masses.array();
      gradientMatrix.col(2).array() /= masses.array();
      if (transformation) {
        gradients = transformation->gradientsToInternal(gradientMatrix);
      }
      else {
        gradients = Eigen::Map<const Eigen::VectorXd>(gradientMatrix.data(), nAtoms * 3);
      }
    };
    // Get initial positions as (transformed) array
    Eigen::VectorXd positions;
    if (transformation) {
      positions = transformation->coordinatesToInternal(coordinates);
    }
    else {
      positions = Eigen::Map<const Eigen::VectorXd>(coordinates.data(), nAtoms * 3);
    }
    // Optimize
    int cycles = 0;
    try {
      cycles = optimizer.optimize(positions, update, check);
    }
    catch (const InternalCoordinatesException& e) {
      log.output << "Internal coordinates broke down. Continuing in Cartesians." << Core::Log::nl;
      // Update coordinates to the last ones that were successfully reconverted
      Utils::PositionCollection lastCoordinates = _calculator.getPositions();
      // Disable coordinate transformation
      this->coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
      transformation = std::make_shared<InternalCoordinates>(atoms, true);
      Eigen::VectorXd lastPositions = transformation->coordinatesToInternal(lastCoordinates);
      // Restart optimization
      optimizer.prepareRestart(cycle);
      cycles = optimizer.optimize(lastPositions, update, check);
      positions = lastPositions;
    }
    // Update AtomCollection and return
    if (transformation) {
      coordinates = transformation->coordinatesToCartesian(positions);
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
  void setSettings(const Settings& settings) override {
    check.applySettings(settings);
    optimizer.applySettings(settings);
    this->initialStepSize = settings.getDouble(IrcOptimizerBase::ircInitialStepSizeKey);
    this->coordinateSystem = CoordinateSystemInterpreter::getCoordinateSystemFromString(
        settings.getString(IrcOptimizerBase::ircCoordinateSystemKey));
  };
  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  Settings getSettings() const override {
    return IRCOptimizerSettings<OptimizerType, GradientBasedCheck>(*this, optimizer, check);
  };
  /**
   * @brief Get the settings of the calculator used for the energy calculations during the optimization.
   * @return std::shared_ptr<Settings> The settings of the calculator.
   */
  const std::shared_ptr<Settings> getCalculatorSettings() const override {
    return std::make_shared<Settings>(_calculator.settings());
  };
  /**
   * @brief Add an observer function that will be triggered in each iteration.
   *
   * @param function A function to be executed in every loop of the optimization.
   *                 The function will have access to the current cycle count
   *                 the current value and to a const reference of the current
   *                 parameters.
   */
  void addObserver(std::function<void(const int&, const double&, const Eigen::VectorXd&)> function) final {
    optimizer.addObserver(function);
  }
  /**
   * @brief Clear all existing observer functions.
   *
   * For optimization problems with very fast evaluations of the underlying function
   * the removal of all observers can increase performance as the observers are given
   * as std::functions and can not be added via templates.
   */
  void clearObservers() final {
    optimizer.clearObservers();
  }
  /**
   * @brief The underlying convergence check
   *
   * @note getter to be accessible via base class
   * @return GradientBasedCheck the class holding all convergence thresholds.
   */
  const GradientBasedCheck getConvergenceCheck() const override {
    return check;
  };
  /// @brief The underlying optimizer, public in order to change it's settings.
  OptimizerType optimizer;
  /// @brief The underlying convergence check, public in order to change it's settings.
  GradientBasedCheck check;

 private:
  Core::Calculator& _calculator;
};

/*=============================*
 *       Gradient Based
 *=============================*/

template<>
inline IrcOptimizer<SteepestDescent>::IrcOptimizer(Core::Calculator& calculator) : _calculator(calculator) {
  // Set SD factor to 2.0
  this->optimizer.factor = 2.0;
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_GEOMETRYOPTIMIZER_H_
