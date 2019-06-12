/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_AFIROPTIMIZER_H_
#define UTILS_AFIROPTIMIZER_H_

#include "Utils/CalculatorBasics/PropertyList.h"
#include "Utils/CalculatorBasics/Results.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Geometry/GeometryUtilities.h"
#include "Utils/GeometryOptimization/AFIROptimizerBase.h"
#include "Utils/GeometryOptimization/AFIROptimizerSettings.h"
#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include <Core/Interfaces/Calculator.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @brief A version of the GeometryOptimizer that optimizes the underlying structure while applying an
 *        additional artificial force for the induction of reactions.
 *
 * For more details about AFIR see the base class Utils::AFIROptimizerBase .
 *
 * @tparam OptimizerType Expects any of the Optimizer classes. Note that some special optimizers
 *                       may not yet be supported or may need additional specialization.
 */
template<class OptimizerType>
class AFIROptimizer : public AFIROptimizerBase {
 public:
  /**
   * @brief Construct a new AFIROptimizer object.
   * @param calculator The calculator to be used for the underlying single point/gradient calculations.
   */
  AFIROptimizer(Core::Calculator& calculator) : _calculator(calculator){};

  /**
   * @brief See AFIROptimizerBase::optimize().
   *
   * @param atoms The AtomCollection (Geometry) to be optimized.
   * @return int  The final number of optimization cycles carried out.
   */
  virtual int optimize(AtomCollection& atoms) final {
    // Configure Calculator
    _calculator.setStructure(atoms);
    _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
    // Transformation into internal basis
    Eigen::MatrixXd transformation;
    if (this->transformCoordinates) {
      transformation = Geometry::calculateRotTransFreeTransformMatrix(atoms.getPositions(), atoms.getElements());
    }
    // Define update function
    int cycle = 0;
    const unsigned int nAtoms = atoms.size();
    auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
      cycle++;
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
      // Evaluate AF
      auto currentStructure = _calculator.getStructure();
      double afEnergy;
      Eigen::MatrixXd afGradients(gradientMatrix);
      this->evaluateArtificialForces(*currentStructure, afEnergy, afGradients);
      // Add AF contributions
      value += afEnergy;
      gradientMatrix += std::min(cycle / (phaseIn + 1.0), 1.0) * afGradients;
      // Store gradients
      if (this->transformCoordinates) {
        auto tmp = Eigen::Map<const Eigen::VectorXd>(gradientMatrix.data(), nAtoms * 3);
        gradients = (transformation.transpose() * tmp).eval();
      }
      else {
        gradients = Eigen::Map<const Eigen::VectorXd>(gradientMatrix.data(), nAtoms * 3);
      }
    };
    // Get initial positions
    Eigen::VectorXd positions;
    if (this->transformCoordinates) {
      auto tmp = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), nAtoms * 3);
      positions = (transformation.transpose() * tmp).eval();
    }
    else {
      positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), nAtoms * 3);
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
    this->rhsList = settings.getIntList(AFIROptimizerBase::afirRHSListKey);
    this->lhsList = settings.getIntList(AFIROptimizerBase::afirLHSListKey);
    this->weak = settings.getBool(AFIROptimizerBase::afirWeakForcesKey);
    this->attractive = settings.getBool(AFIROptimizerBase::afirAttractiveKey);
    this->energyAllowance = settings.getDouble(AFIROptimizerBase::afirEnergyAllowanceKey);
    this->phaseIn = settings.getInt(AFIROptimizerBase::afirPhaseInKey);
    this->transformCoordinates = settings.getBool(AFIROptimizerBase::afirTransfromCoordinatesKey);
  };
  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  virtual Settings getSettings() const override {
    return AFIROptimizerSettings<OptimizerType, GradientBasedCheck>(*this, optimizer, check);
  };
  /**
   * @brief Add an observer function that will be triggered in each iteration.
   *
   * @param function A function to be executed in every loop of the optimization.
   *                 The function will have access to the current cycle count,
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

#endif // UTILS_AFIROPTIMIZER_H_
