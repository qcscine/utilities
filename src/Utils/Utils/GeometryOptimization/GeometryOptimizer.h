/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_GEOMETRYOPTIMIZER_H_
#define UTILS_GEOMETRYOPTIMIZER_H_

#include "CoordinateSystem.h"
#include "Utils/CalculatorBasics.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/InternalCoordinates.h"
#include "Utils/IO/ChemicalFileFormats/XyzStreamHandler.h"
#include "Utils/Optimizer/GradientBased/Bfgs.h"
#include "Utils/Optimizer/GradientBased/Dimer.h"
#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
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
 * @brief The base class for all Geometry optimizers.
 *
 * The purpose of the geometry optimizers is to wrap technical details needed
 * for the actual optimization into a class that is easily used and has a small
 * set of data needed in its constructor, is mainly configured via its settings,
 * and then exposes the optimization through a simple function accepting a geometry
 * (AtomCollection)\n
 * \n
 * The main purpose of this base class is to hide the template parameter(s)
 * of the derived class(es).
 */
class GeometryOptimizerBase {
 public:
  static constexpr const char* geooptFixedAtomsKey = "geoopt_constrained_atoms";
  static constexpr const char* geooptCoordinateSystemKey = "geoopt_coordinate_system";
  /// @brief Default constructor.
  GeometryOptimizerBase() = default;
  /// @brief Virtual default destructor.
  virtual ~GeometryOptimizerBase() = default;
  /**
   * @brief The main functionality of the geometry optimizer.
   *
   * This function wraps the optimize functions of the underlying optimizer.
   *
   * @param atoms The AtomCollection (Geometry) to be optimized.
   * @param log The logger to which eventual output is written.
   * @return int  The final number of optimization cycles carried out.
   */
  virtual int optimize(AtomCollection& atoms, Core::Log& log) = 0;
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
  /**
   * @brief Get a copy of the settings of the calculator used for the energy calculations during the optimization.
   * @return std::shared_ptr<Settings> The settings of the calculator.
   */
  virtual std::shared_ptr<Settings> getCalculatorSettings() const = 0;
  /**
   * @brief The underlying convergence check
   *
   * @return GradientBasedCheck the class holding all convergence thresholds.
   */
  virtual const GradientBasedCheck& getConvergenceCheck() const = 0;
  /**
   * @brief Vector containing the atom indices to which Cartesian constraints are applied during the optimization.
   *        Note that only an empty vector still allows for the use of internal coordinates.
   */
  std::vector<int> fixedAtoms = {};

  /**
   * @brief Set the coordinate system in which the optimization shall be performed
   *
   * The optimization can be carried out in the internal coordinate space or with removed translations and rotations
   * possibly accelerating convergence.
   */
  CoordinateSystem coordinateSystem = CoordinateSystem::Internal;
};

/**
 * @brief Settings for a GeometryOptimizer
 *
 * Uses template arguments in order to automatically include the
 * settings of underlying objects into the given settings.
 *
 * @tparam OptimizerType The underlying Optimizer class.
 * @tparam ConvergenceCheckType The underlying ConvergenceCheck class.
 */
template<class OptimizerType, class ConvergenceCheckType>
class GeometryOptimizerSettings : public Settings {
 public:
  /**
   * @brief Construct a new GeometryOptimizerSettings object.
   *
   * Sets the default values of the settings to the current values set in the objects
   * given to the constructor.
   *
   * @param base The geometry optimizer.
   * @param optimizer The optimizer.
   * @param check The convergence check criteria.
   */
  GeometryOptimizerSettings(const GeometryOptimizerBase& base, const OptimizerType& optimizer, const ConvergenceCheckType& check)
    : Settings("GeometryOptimizerSettings") {
    optimizer.addSettingsDescriptors(this->_fields);
    check.addSettingsDescriptors(this->_fields);

    UniversalSettings::OptionListDescriptor geoopt_coordinate_system("Set the coordinate system.");
    geoopt_coordinate_system.addOption("internal");
    geoopt_coordinate_system.addOption("cartesianWithoutRotTrans");
    geoopt_coordinate_system.addOption("cartesian");
    geoopt_coordinate_system.setDefaultOption(CoordinateSystemInterpreter::getStringFromCoordinateSystem(base.coordinateSystem));
    this->_fields.push_back(GeometryOptimizerBase::geooptCoordinateSystemKey, std::move(geoopt_coordinate_system));

    UniversalSettings::IntListDescriptor geooptFixedAtoms(
        "List of atoms with Cartesian constraints applied to them during the optimization.");
    this->_fields.push_back(GeometryOptimizerBase::geooptFixedAtomsKey, std::move(geooptFixedAtoms));

    this->resetToDefaults();
  }
};

/**
 * @brief Basically just a templated version of the base class GeometryOptimizerBase,
 *        where the template defines the actual optimizer used in the geometry optimization.
 *
 * @tparam OptimizerType Expects any of the Optimizer classes. Note that some special optimizers
 *                       may not yet be supported or may need additional specialization.
 */
template<class OptimizerType>
class GeometryOptimizer : public GeometryOptimizerBase {
 public:
  /**
   * @brief Construct a new GeometryOptimizer object.
   * @param calculator The calculator to be used for the single point/gradient calculations.
   */
  GeometryOptimizer(Core::Calculator& calculator) : _calculator(calculator) {
    // With current internal coordinates not possible
    if (std::is_same<OptimizerType, Lbfgs>::value || std::is_same<OptimizerType, NewtonRaphson>::value ||
        std::is_same<OptimizerType, Bofill>::value || std::is_same<OptimizerType, EigenvectorFollowing>::value) {
      this->coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
    }
    // Performs better with Cartesian coordinates with current internal coordinates
    if (std::is_same<OptimizerType, Dimer>::value) {
      this->coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
    }
  };
  /**
   * @brief See GeometryOptimizerBase::optimize().
   * @param log The logger to which eventual output is written.
   * @param atoms The AtomCollection (Geometry) to be optimized.
   *
   * @return int  The final number of optimization cycles carried out.
   */
  int optimize(AtomCollection& atoms, Core::Log& log) final {
    // Disable L-BFGS + internals for now
    // TODO: fix the Hessian projection in the L-BFGS to allow for this combination
    if (std::is_same<OptimizerType, Lbfgs>::value && this->coordinateSystem == CoordinateSystem::Internal) {
      throw std::runtime_error("Error: The L-BFGS optimizer is currently not allowed with Internal coordinates.");
    }
    // Configure Calculator
    _calculator.setStructure(atoms);
    _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
    // Transformation into internal coordinates
    std::shared_ptr<InternalCoordinates> transformation = nullptr;
    if (this->coordinateSystem == CoordinateSystem::Internal) {
      transformation = std::make_shared<InternalCoordinates>(atoms);
    }
    else if (this->coordinateSystem == CoordinateSystem::CartesianWithoutRotTrans) {
      transformation = std::make_shared<InternalCoordinates>(atoms, true);
    }
    // Define update function
    const unsigned int nAtoms = atoms.size();
    // Count cycles for eventual restart
    int cycle = 0;
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
      _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
      atoms.setPositions(coordinates);
      Results results =
          CalculationRoutines::calculateWithCatch(_calculator, log, "Aborting optimization due to failed calculation");
      value = results.get<Property::Energy>();
      if (transformation) {
        gradients = transformation->gradientsToInternal(results.get<Property::Gradients>());
      }
      else {
        gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
      }
    };
    // Get initial positions
    Eigen::VectorXd positions;
    if (transformation) {
      positions = transformation->coordinatesToInternal(atoms.getPositions());
    }
    else {
      positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), atoms.size() * 3);
      if (!fixedAtoms.empty()) {
        optimizer.mask.resizeLike(positions);
        optimizer.mask.setConstant(true);
        for (const int fixedAtom : fixedAtoms) {
          if (fixedAtom < 0 || fixedAtom >= static_cast<int>(nAtoms)) {
            throw std::runtime_error("Constrained atom index " + std::to_string(fixedAtom) + " is invalid!");
          }
          optimizer.mask.template segment<3>(3 * fixedAtom).setConstant(false);
        }
      }
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
      // Disable true internal coordinates
      this->coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
      transformation = std::make_shared<InternalCoordinates>(atoms, true);
      Eigen::VectorXd lastPositions = transformation->coordinatesToInternal(lastCoordinates);
      // Restart optimization
      optimizer.prepareRestart(cycle);
      cycles = optimizer.optimize(lastPositions, update, check);
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
    optimizer.applySettings(settings);
    this->coordinateSystem = CoordinateSystemInterpreter::getCoordinateSystemFromString(
        settings.getString(GeometryOptimizerBase::geooptCoordinateSystemKey));

    // For Cartesian constraints:
    this->fixedAtoms = settings.getIntList(GeometryOptimizerBase::geooptFixedAtomsKey);

    // Check whether constraints and coordinate transformations are both switched on:
    if (!this->fixedAtoms.empty() && this->coordinateSystem != CoordinateSystem::Cartesian) {
      throw std::logic_error("Cartesian constraints cannot be set when using coordinate transformations! Set "
                             "'geoopt_coordinate_system' to 'cartesian'.");
    }
  };
  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  Settings getSettings() const override {
    return GeometryOptimizerSettings<OptimizerType, GradientBasedCheck>(*this, optimizer, check);
  };
  /**
   * @brief Get a copy of the settings of the calculator used for the energy calculations during the optimization.
   * @return std::shared_ptr<Settings> The settings of the calculator.
   */
  std::shared_ptr<Settings> getCalculatorSettings() const override {
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
  const GradientBasedCheck& getConvergenceCheck() const override {
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
 *  Approximate Hessian Based
 *=============================*/

template<>
inline int GeometryOptimizer<Dimer>::optimize(AtomCollection& atoms, Core::Log& log) {
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
  // Define update function
  const unsigned int nAtoms = atoms.size();
  auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
    Utils::PositionCollection coordinates;
    if (transformation) {
      coordinates = transformation->coordinatesToCartesian(parameters);
    }
    else {
      coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
    }
    _calculator.modifyPositions(coordinates);
    _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
    atoms.setPositions(coordinates);
    Results results =
        CalculationRoutines::calculateWithCatch(_calculator, log, "Aborting optimization due to failed calculation");
    value = results.get<Property::Energy>();
    if (transformation) {
      gradients = transformation->gradientsToInternal(results.get<Property::Gradients>());
    }
    else {
      gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
    }
  };
  // Get initial positions
  Eigen::VectorXd positions;
  if (transformation) {
    positions = transformation->coordinatesToInternal(atoms.getPositions());
  }
  else {
    positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), atoms.size() * 3);
    if (!fixedAtoms.empty()) {
      optimizer.mask.resizeLike(positions);
      optimizer.mask.setConstant(true);
      for (const int fixedAtom : fixedAtoms) {
        if (fixedAtom < 0 || fixedAtom >= static_cast<int>(nAtoms)) {
          throw std::runtime_error("Constrained atom index " + std::to_string(fixedAtom) + " is invalid!");
        }
        optimizer.mask.template segment<3>(3 * fixedAtom).setConstant(false);
      }
    }
  }
  // Get projection
  if (this->coordinateSystem == CoordinateSystem::Internal) {
    optimizer.invH = transformation->inverseHessianGuess();
    optimizer.projection = std::make_shared<std::function<void(Eigen::MatrixXd&)>>(
        [&transformation](Eigen::MatrixXd& inv) { inv = transformation->projectHessianInverse(inv); });
  }
  else {
    optimizer.projection = nullptr;
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
    // Disable true internal coordinates
    this->coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
    transformation = std::make_shared<InternalCoordinates>(atoms, true);
    Eigen::VectorXd lastPositions = transformation->coordinatesToInternal(lastCoordinates);
    // Restart optimization
    // Get cycle from optimizer because one cycle contains several calls to the update function
    optimizer.prepareRestart(optimizer.getCycle());
    cycles = optimizer.optimize(lastPositions, update, check);
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

template<>
inline int GeometryOptimizer<Bfgs>::optimize(AtomCollection& atoms, Core::Log& log) {
  // Reset Bfgs
  this->optimizer.invH.resize(0, 0);
  // Configure Calculator
  _calculator.setStructure(atoms);
  _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  // Transformation into internal coordinates
  std::shared_ptr<InternalCoordinates> transformation = nullptr;
  if (this->coordinateSystem == CoordinateSystem::Internal) {
    transformation = std::make_shared<InternalCoordinates>(atoms);
  }
  else if (this->coordinateSystem == CoordinateSystem::CartesianWithoutRotTrans) {
    transformation = std::make_shared<InternalCoordinates>(atoms, true);
  }
  // Define update function
  const unsigned int nAtoms = atoms.size();
  // Count cycles for eventual restart
  int cycle = 0;
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
    _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
    atoms.setPositions(coordinates);
    Results results =
        CalculationRoutines::calculateWithCatch(_calculator, log, "Aborting optimization due to failed calculation");
    value = results.get<Property::Energy>();
    if (transformation) {
      gradients = transformation->gradientsToInternal(results.get<Property::Gradients>());
    }
    else {
      gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
    }
  };
  // Get initial positions
  Eigen::VectorXd positions;
  if (transformation) {
    positions = transformation->coordinatesToInternal(atoms.getPositions());
  }
  else {
    positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), atoms.size() * 3);
    if (!fixedAtoms.empty()) {
      optimizer.mask.resizeLike(positions);
      optimizer.mask.setConstant(true);
      for (const int fixedAtom : fixedAtoms) {
        if (fixedAtom < 0 || fixedAtom >= static_cast<int>(nAtoms)) {
          throw std::runtime_error("Constrained atom index " + std::to_string(fixedAtom) + " is invalid!");
        }
        optimizer.mask.template segment<3>(3 * fixedAtom).setConstant(false);
      }
    }
  }
  // Get projection
  if (this->coordinateSystem == CoordinateSystem::Internal) {
    optimizer.invH = transformation->inverseHessianGuess();
    optimizer.projection = [&transformation](Eigen::MatrixXd& inv) { inv = transformation->projectHessianInverse(inv); };
  }
  else {
    optimizer.projection = nullptr;
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
    // Disable true internal coordinates
    this->coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
    transformation = std::make_shared<InternalCoordinates>(atoms, true);
    Eigen::VectorXd lastPositions = transformation->coordinatesToInternal(lastCoordinates);
    // Restart optimization
    optimizer.prepareRestart(cycle);
    cycles = optimizer.optimize(lastPositions, update, check);
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

/*================================*
 *  Hessian Base Specializations
 *================================*/
template<>
inline int GeometryOptimizer<Bofill>::optimize(AtomCollection& atoms, Core::Log& log) {
  if (this->coordinateSystem == CoordinateSystem::Internal) {
    throw std::runtime_error("Error: The Bofill optimizer is currently not allowed with Internal coordinates.");
  }
  // Configure Calculator
  _calculator.setStructure(atoms);
  _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
  // Transformation into internal basis
  std::shared_ptr<InternalCoordinates> transformation = nullptr;
  if (this->coordinateSystem == CoordinateSystem::CartesianWithoutRotTrans) {
    transformation = std::make_shared<InternalCoordinates>(atoms, true);
  }
  // Define update function
  const unsigned int nAtoms = atoms.size();
  auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                          Eigen::MatrixXd& hessian, bool calcHessian) {
    Utils::PositionCollection coordinates;
    if (transformation) {
      coordinates = transformation->coordinatesToCartesian(parameters);
    }
    else {
      coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
    }
    _calculator.modifyPositions(coordinates);
    atoms.setPositions(coordinates);

    if (calcHessian) {
      _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
      Results results =
          CalculationRoutines::calculateWithCatch(_calculator, log, "Aborting optimization due to failed calculation");
      value = results.get<Property::Energy>();
      if (transformation) {
        gradients = transformation->gradientsToInternal(results.get<Property::Gradients>());
        hessian = transformation->hessianToInternal(results.get<Property::Hessian>());
      }
      else {
        gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
        hessian = results.get<Property::Hessian>();
      }
    }
    else {
      _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
      Results results =
          CalculationRoutines::calculateWithCatch(_calculator, log, "Aborting optimization due to failed calculation");
      value = results.get<Property::Energy>();
      if (transformation) {
        gradients = transformation->gradientsToInternal(results.get<Property::Gradients>());
      }
      else {
        gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
      }
    }
  };
  // Get initial positions
  Eigen::VectorXd positions;
  if (transformation) {
    positions = transformation->coordinatesToInternal(atoms.getPositions());
  }
  else {
    positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), nAtoms * 3);
    if (!fixedAtoms.empty()) {
      optimizer.mask.resizeLike(positions);
      optimizer.mask.setConstant(true);
      for (const int fixedAtom : fixedAtoms) {
        if (fixedAtom < 0 || fixedAtom >= static_cast<int>(nAtoms)) {
          throw std::runtime_error("Constrained atom index " + std::to_string(fixedAtom) + " is invalid!");
        }
        optimizer.mask.template segment<3>(3 * fixedAtom).setConstant(false);
      }
    }
  }

  // Optimize
  auto cycles = optimizer.optimize(positions, update, check);
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

template<>
inline int GeometryOptimizer<NewtonRaphson>::optimize(AtomCollection& atoms, Core::Log& log) {
  if (this->coordinateSystem == CoordinateSystem::Internal) {
    throw std::runtime_error("Error: The Newton Raphson optimizer is currently not allowed with Internal coordinates.");
  }
  // Configure Calculator
  _calculator.setStructure(atoms);
  _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
  // Transformation into internal basis
  std::shared_ptr<InternalCoordinates> transformation = nullptr;
  if (this->coordinateSystem == CoordinateSystem::CartesianWithoutRotTrans) {
    transformation = std::make_shared<InternalCoordinates>(atoms, true);
  }
  // Define update function
  const unsigned int nAtoms = atoms.size();
  auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                          Eigen::MatrixXd& hessian, bool /*calcHessian*/) {
    Utils::PositionCollection coordinates;
    if (transformation) {
      coordinates = transformation->coordinatesToCartesian(parameters);
    }
    else {
      coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
    }
    _calculator.modifyPositions(coordinates);
    atoms.setPositions(coordinates);
    _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
    Results results =
        CalculationRoutines::calculateWithCatch(_calculator, log, "Aborting optimization due to failed calculation");
    value = results.get<Property::Energy>();
    if (transformation) {
      gradients = transformation->gradientsToInternal(results.get<Property::Gradients>());
      hessian = transformation->hessianToInternal(results.get<Property::Hessian>());
    }
    else {
      gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
      hessian = results.get<Property::Hessian>();
    }
  };
  // Get initial positions
  Eigen::VectorXd positions;
  if (transformation) {
    positions = transformation->coordinatesToInternal(atoms.getPositions());
  }
  else {
    positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), nAtoms * 3);
    if (!fixedAtoms.empty()) {
      optimizer.mask.resizeLike(positions);
      optimizer.mask.setConstant(true);
      for (const int fixedAtom : fixedAtoms) {
        if (fixedAtom < 0 || fixedAtom >= static_cast<int>(nAtoms)) {
          throw std::runtime_error("Constrained atom index " + std::to_string(fixedAtom) + " is invalid!");
        }
        optimizer.mask.template segment<3>(3 * fixedAtom).setConstant(false);
      }
    }
  }

  // Optimize
  auto cycles = optimizer.optimize(positions, update, check);
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

template<>
inline int GeometryOptimizer<EigenvectorFollowing>::optimize(AtomCollection& atoms, Core::Log& log) {
  if (this->coordinateSystem == CoordinateSystem::Internal) {
    throw std::runtime_error(
        "Error: The EigenvectorFollowing optimizer is currently not allowed with Internal coordinates.");
  }
  // Configure Calculator
  _calculator.setStructure(atoms);
  _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
  // Transformation into internal basis
  std::shared_ptr<InternalCoordinates> transformation = nullptr;
  if (this->coordinateSystem == CoordinateSystem::CartesianWithoutRotTrans) {
    transformation = std::make_shared<InternalCoordinates>(atoms, true);
  }
  // Define update function
  const unsigned int nAtoms = atoms.size();
  auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                          Eigen::MatrixXd& hessian, bool /*calcHessian*/) {
    Utils::PositionCollection coordinates;
    if (transformation) {
      coordinates = transformation->coordinatesToCartesian(parameters);
    }
    else {
      coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
    }
    _calculator.modifyPositions(coordinates);
    atoms.setPositions(coordinates);
    _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
    Results results =
        CalculationRoutines::calculateWithCatch(_calculator, log, "Aborting optimization due to failed calculation");
    value = results.get<Property::Energy>();
    if (transformation) {
      gradients = transformation->gradientsToInternal(results.get<Property::Gradients>());
      hessian = transformation->hessianToInternal(results.get<Property::Hessian>());
    }
    else {
      gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
      hessian = results.get<Property::Hessian>();
    }
  };
  // Get initial positions
  Eigen::VectorXd positions;
  if (transformation) {
    positions = transformation->coordinatesToInternal(atoms.getPositions());
  }
  else {
    positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), nAtoms * 3);
    if (!fixedAtoms.empty()) {
      optimizer.mask.resizeLike(positions);
      optimizer.mask.setConstant(true);
      for (const int fixedAtom : fixedAtoms) {
        if (fixedAtom < 0 || fixedAtom >= static_cast<int>(nAtoms)) {
          throw std::runtime_error("Constrained atom index " + std::to_string(fixedAtom) + " is invalid!");
        }
        optimizer.mask.template segment<3>(3 * fixedAtom).setConstant(false);
      }
    }
  }

  // Optimize
  auto cycles = optimizer.optimize(positions, update, check);
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

#endif // UTILS_GEOMETRYOPTIMIZER_H_
