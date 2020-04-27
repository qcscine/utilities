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
#include "Utils/Geometry/InternalCoordinates.h"
#include "Utils/GeometryOptimization/AfirOptimizerBase.h"
#include "Utils/GeometryOptimization/AfirOptimizerSettings.h"
#include "Utils/Optimizer/GradientBased/Bfgs.h"
#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/GradientBased/Lbfgs.h"
#include <Core/Interfaces/Calculator.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @brief A version of the GeometryOptimizer that optimizes the underlying structure while applying an
 *        additional artificial force for the induction of reactions.
 *
 * For more details about AFIR see the base class Utils::AfirOptimizerBase .
 *
 * @tparam OptimizerType Expects any of the Optimizer classes. Note that some special optimizers
 *                       may not yet be supported or may need additional specialization.
 */
template<class OptimizerType>
class AfirOptimizer : public AfirOptimizerBase {
 public:
  /**
   * @brief Construct a new AfirOptimizer object.
   * @param calculator The calculator to be used for the underlying single point/gradient calculations.
   */
  AfirOptimizer(Core::Calculator& calculator) : _calculator(calculator){};

  /**
   * @brief See AfirOptimizerBase::optimize().
   *
   * @param atoms The AtomCollection (Geometry) to be optimized.
   * @return int  The final number of optimization cycles carried out.
   */
  virtual int optimize(AtomCollection& atoms) final {
    // Disable L-BFGS + internals for now
    //  TODO fix the hessian projection in the L-BFGS to allow for this combination
    if (std::is_same<OptimizerType, Lbfgs>::value && this->transformCoordinates)
      throw std::runtime_error("Error: L-BFGS + Internal coordinates are currently not allowed.");
    // Configure Calculator
    _calculator.setStructure(atoms);
    _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
    // Transform into internal coordinates
    auto transformation = std::make_shared<InternalCoordinates>(atoms);
    // Define update function
    int cycle = 0;
    const unsigned int nAtoms = atoms.size();
    auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
      cycle++;
      Utils::PositionCollection coordinates;
      if (this->transformCoordinates) {
        coordinates = transformation->coordinatesToCartesian(parameters);
      }
      else {
        coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
      }
      _calculator.modifyPositions(coordinates);
      Utils::Results results = _calculator.calculate("Geometry Optimization Cycle");
      value = results.get<Property::Energy>();
      auto gradientMatrix = results.get<Property::Gradients>();
      // Evaluate AF
      auto currentStructure = _calculator.getStructure();
      double afEnergy;
      Utils::GradientCollection afGradients(gradientMatrix);
      this->evaluateArtificialForces(*currentStructure, afEnergy, afGradients);
      // Add AF contributions
      value += afEnergy;
      gradientMatrix += std::min(cycle / (phaseIn + 1.0), 1.0) * afGradients;
      // Store gradients
      if (this->transformCoordinates) {
        gradients = transformation->gradientsToInternal(gradientMatrix);
      }
      else {
        gradients = Eigen::Map<const Eigen::VectorXd>(gradientMatrix.data(), nAtoms * 3);
      }
    };
    // Get initial positions
    Eigen::VectorXd positions;
    if (this->transformCoordinates) {
      positions = transformation->coordinatesToInternal(atoms.getPositions());
    }
    else {
      positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), atoms.size() * 3);
    }
    // Optimize
    auto cycles = optimizer.optimize(positions, update, check);
    // Update Atom collection and return
    Utils::PositionCollection coordinates;
    if (this->transformCoordinates) {
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
  virtual void setSettings(const Settings& settings) override {
    check.applySettings(settings);
    optimizer.applySettings(settings);
    this->rhsList = settings.getIntList(AfirOptimizerBase::afirRHSListKey);
    this->lhsList = settings.getIntList(AfirOptimizerBase::afirLHSListKey);
    this->weak = settings.getBool(AfirOptimizerBase::afirWeakForcesKey);
    this->attractive = settings.getBool(AfirOptimizerBase::afirAttractiveKey);
    this->energyAllowance = settings.getDouble(AfirOptimizerBase::afirEnergyAllowanceKey);
    this->phaseIn = settings.getInt(AfirOptimizerBase::afirPhaseInKey);
    this->transformCoordinates = settings.getBool(AfirOptimizerBase::afirTransfromCoordinatesKey);
  };
  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  virtual Settings getSettings() const override {
    return AfirOptimizerSettings<OptimizerType, GradientBasedCheck>(*this, optimizer, check);
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

/*===================*
 *  Specializations
 *===================*/

template<>
inline int AfirOptimizer<Bfgs>::optimize(AtomCollection& atoms) {
  // Configure Calculator
  _calculator.setStructure(atoms);
  _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  // Transform into internal coordinates
  auto transformation = std::make_shared<InternalCoordinates>(atoms);
  // Define update function
  int cycle = 0;
  const unsigned int nAtoms = atoms.size();
  auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
    cycle++;
    Utils::PositionCollection coordinates;
    if (this->transformCoordinates) {
      coordinates = transformation->coordinatesToCartesian(parameters);
    }
    else {
      coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
    }
    _calculator.modifyPositions(coordinates);
    Utils::Results results = _calculator.calculate("Geometry Optimization Cycle");
    value = results.get<Property::Energy>();
    auto gradientMatrix = results.get<Property::Gradients>();
    // Evaluate AF
    auto currentStructure = _calculator.getStructure();
    double afEnergy;
    Utils::GradientCollection afGradients(gradientMatrix);
    this->evaluateArtificialForces(*currentStructure, afEnergy, afGradients);
    // Add AF contributions
    value += afEnergy;
    gradientMatrix += std::min(cycle / (phaseIn + 1.0), 1.0) * afGradients;
    // Store gradients
    if (this->transformCoordinates) {
      gradients = transformation->gradientsToInternal(gradientMatrix);
    }
    else {
      gradients = Eigen::Map<const Eigen::VectorXd>(gradientMatrix.data(), nAtoms * 3);
    }
  };
  // Get initial positions
  Eigen::VectorXd positions;
  if (this->transformCoordinates) {
    positions = transformation->coordinatesToInternal(atoms.getPositions());
  }
  else {
    positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), atoms.size() * 3);
  }
  // Optimize
  if (this->transformCoordinates) {
    optimizer.projection = std::make_unique<std::function<void(Eigen::MatrixXd&)>>(
        [&transformation](Eigen::MatrixXd& inv) { inv = transformation->projectHessianInverse(inv); });
  }
  auto cycles = optimizer.optimize(positions, update, check);
  // Update Atom collection and return
  Utils::PositionCollection coordinates;
  if (this->transformCoordinates) {
    coordinates = transformation->coordinatesToCartesian(positions);
  }
  else {
    coordinates = Eigen::Map<const Utils::PositionCollection>(positions.data(), nAtoms, 3);
  }
  atoms.setPositions(coordinates);
  return cycles;
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_AFIROPTIMIZER_H_
