/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_AFIROPTIMIZER_H_
#define UTILS_AFIROPTIMIZER_H_

#include "Utils/CalculatorBasics.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Geometry/GeometryUtilities.h"
#include "Utils/Geometry/InternalCoordinates.h"
#include "Utils/GeometryOptimization/AfirConvergenceCheck.h"
#include "Utils/GeometryOptimization/AfirOptimizerBase.h"
#include "Utils/GeometryOptimization/AfirOptimizerSettings.h"
#include "Utils/Optimizer/GradientBased/Bfgs.h"
#include "Utils/Optimizer/GradientBased/Dimer.h"
#include "Utils/Optimizer/GradientBased/Lbfgs.h"
#include "Utils/Optimizer/HessianBased/Bofill.h"
#include "Utils/Optimizer/HessianBased/EigenvectorFollowing.h"
#include "Utils/Optimizer/HessianBased/NewtonRaphson.h"
#include "Utils/UniversalSettings/SettingsNames.h"
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
   * @param log The logger to which eventual output is written.
   * @return int  The final number of optimization cycles carried out.
   */
  int optimize(AtomCollection& atoms, Core::Log& log) final {
    // Ensure that AFIR check uses the same LHS and RHS list and knows whether internal coordinates are used
    check.lhsList = this->lhsList;
    check.rhsList = this->rhsList;
    // Disable L-BFGS + internals for now
    // TODO fix the Hessian projection in the L-BFGS to allow for this combination
    if (std::is_same<OptimizerType, Lbfgs>::value && this->coordinateSystem == CoordinateSystem::Internal) {
      throw std::runtime_error("Error: This optimizer is currently not allowed with Internal coordinates.");
    }
    if (std::is_same<OptimizerType, NewtonRaphson>::value || std::is_same<OptimizerType, Bofill>::value ||
        std::is_same<OptimizerType, EigenvectorFollowing>::value || std::is_same<OptimizerType, Dimer>::value) {
      throw std::runtime_error("Error: This optimizer is currently not available for AFIR optimization.");
    }
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
    check.transformation = transformation;
    // Define update function
    int cycle = 0;
    const unsigned int nAtoms = atoms.size();
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
      // Evaluate AF
      auto currentStructure = _calculator.getStructure();
      double afEnergy;
      Utils::GradientCollection afGradients(gradientMatrix);
      this->evaluateArtificialForces(*currentStructure, afEnergy, afGradients);
      // Add AF contributions
      value += afEnergy;
      gradientMatrix += std::min(cycle / (phaseIn + 1.0), 1.0) * afGradients;
      // Store gradients
      if (transformation) {
        gradients = transformation->gradientsToInternal(gradientMatrix);
      }
      else {
        gradients = Eigen::Map<const Eigen::VectorXd>(gradientMatrix.data(), nAtoms * 3);
      }
    };
    // Get initial positions
    Eigen::VectorXd positions;
    if (transformation) {
      positions = transformation->coordinatesToInternal(atoms.getPositions());
    }
    else {
      positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), atoms.size() * 3);
    }
    // Optimize
    int cycles = 0;
    try {
      cycles = optimizer.optimize(positions, update, check, log);
    }
    catch (const InternalCoordinatesException& e) {
      log.output << "Internal coordinates broke down. Continuing in Cartesians." << Core::Log::nl;
      // Update coordinates to the last ones that were successfully reconverted
      Utils::PositionCollection lastCoordinates = _calculator.getPositions();
      // Disable coordinate transformation
      this->coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
      transformation = std::make_shared<InternalCoordinates>(atoms, true);
      Eigen::VectorXd lastPositions = transformation->coordinatesToInternal(lastCoordinates);
      check.transformation = transformation;
      // Restart optimization
      optimizer.prepareRestart(cycle);
      cycles = optimizer.optimize(lastPositions, update, check, log);
      positions = lastPositions;
    }
    // Update Atom collection and return
    Utils::PositionCollection coordinates;
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
    check.applyAfirSettings(settings);
    optimizer.applySettings(settings);
    this->rhsList = settings.getIntList(SettingsNames::Optimizations::Afir::rHSList);
    this->lhsList = settings.getIntList(SettingsNames::Optimizations::Afir::lHSList);
    this->weak = settings.getBool(SettingsNames::Optimizations::Afir::weakForces);
    this->attractive = settings.getBool(SettingsNames::Optimizations::Afir::attractive);
    this->energyAllowance = settings.getDouble(SettingsNames::Optimizations::Afir::energyAllowance);
    this->phaseIn = settings.getInt(SettingsNames::Optimizations::Afir::phaseIn);
    this->coordinateSystem = CoordinateSystemInterpreter::getCoordinateSystemFromString(
        settings.getString(SettingsNames::Optimizations::Afir::coordinateSystem));
  };
  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  Settings getSettings() const override {
    return AfirOptimizerSettings<OptimizerType, AfirConvergenceCheck>(*this, optimizer, check);
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
   *                 The function will have access to the current cycle count,
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
  AfirConvergenceCheck check;

 private:
  Core::Calculator& _calculator;
};

/*===================*
 *  Specializations
 *===================*/

template<>
inline int AfirOptimizer<Bfgs>::optimize(AtomCollection& atoms, Core::Log& log) {
  // Ensure that AFIR check uses the same LHS and RHS list
  check.lhsList = this->lhsList;
  check.rhsList = this->rhsList;
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
  check.transformation = transformation;
  // Define update function
  int cycle = 0;
  const unsigned int nAtoms = atoms.size();
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
    // Evaluate AF
    auto currentStructure = _calculator.getStructure();
    double afEnergy;
    Utils::GradientCollection afGradients(gradientMatrix);
    this->evaluateArtificialForces(*currentStructure, afEnergy, afGradients);
    // Add AF contributions
    value += afEnergy;
    gradientMatrix += std::min(cycle / (phaseIn + 1.0), 1.0) * afGradients;
    // Store gradients
    if (transformation) {
      gradients = transformation->gradientsToInternal(gradientMatrix);
    }
    else {
      gradients = Eigen::Map<const Eigen::VectorXd>(gradientMatrix.data(), nAtoms * 3);
    }
  };
  // Get initial positions
  Eigen::VectorXd positions;
  if (transformation) {
    positions = transformation->coordinatesToInternal(atoms.getPositions());
  }
  else {
    positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), atoms.size() * 3);
  }
  // Get projection
  if (this->coordinateSystem == CoordinateSystem::Internal) {
    optimizer.projection = [&transformation](Eigen::MatrixXd& inv) { inv = transformation->projectHessianInverse(inv); };
    optimizer.invH = transformation->inverseHessianGuess();
  }
  else {
    optimizer.projection = nullptr;
  }
  // Optimize
  int cycles = 0;
  try {
    cycles = optimizer.optimize(positions, update, check, log);
  }
  catch (const InternalCoordinatesException& e) {
    log.output << "Internal coordinates broke down. Continuing in cartesians." << Core::Log::nl;
    // Update coordinates to the last ones that were successfully reconverted
    Utils::PositionCollection lastCoordinates = _calculator.getPositions();
    // Disable coordinate transformation
    this->coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
    transformation = nullptr;
    check.transformation = nullptr;
    // Restart optimization
    optimizer.prepareRestart(cycle);
    Eigen::VectorXd lastPositions = Eigen::Map<const Eigen::VectorXd>(lastCoordinates.data(), atoms.size() * 3);
    cycles = optimizer.optimize(lastPositions, update, check, log);
    positions = lastPositions;
  }
  // Update Atom collection and return
  Utils::PositionCollection coordinates;
  if (transformation) {
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
