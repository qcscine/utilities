/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_UNITCELLGEOMETRYOPTIMIZER_H_
#define UTILS_UNITCELLGEOMETRYOPTIMIZER_H_

#include "CoordinateSystem.h"
#include "Utils/CalculatorBasics.h"
#include "Utils/DataStructures/PeriodicBoundaries.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/InternalCoordinates.h"
#include "Utils/Geometry/PeriodicSystem.h"
#include "Utils/GeometryOptimization/GeometryOptimizer.h"
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
class UnitCellGeometryOptimizerBase : public GeometryOptimizerBase {
 public:
  static constexpr const char* optimizeAnglesKey = "cellopt_optimize_angles";
  static constexpr const char* optimizeAKey = "cellopt_optimize_a";
  static constexpr const char* optimizeBKey = "cellopt_optimize_b";
  static constexpr const char* optimizeCKey = "cellopt_optimize_c";
  static constexpr const char* geooptMaxIterations = "cellopt_geoopt_max_convergence_iterations";
  static constexpr const char* celloptMaxIterations = "cellopt_cellopt_max_convergence_iterations";

  /// @brief Default constructor.
  UnitCellGeometryOptimizerBase() = default;

  /// @brief Virtual default destructor.
  virtual ~UnitCellGeometryOptimizerBase() = default;

  /**
   * @brief The main functionality of the UnitCellGeometryOptimizer.
   *
   * This function wraps the optimize functions of the underlying optimizers.
   * It has several overloads which take the periodic boundaries from different
   * data structures.
   *
   * @param atoms The AtomCollection to be optimized, the pbc it taken from the calculator
   * and updated in the calculator after the optimization.
   * @param log The logger to which eventual output is written.
   * @return int  The final number of optimization cycles carried out.
   */
  virtual int optimize(AtomCollection& system, Core::Log& log) = 0;
  /**
   * @brief The main functionality of the UnitCellGeometryOptimizer.
   *
   * This function wraps the optimize functions of the underlying optimizers.
   * It has several overloads which take the periodic boundaries from different
   * data structures.
   *
   * @param system The PeriodicSystem to be optimized.
   * @param log The logger to which eventual output is written.
   * @return int  The final number of optimization cycles carried out.
   */
  virtual int optimize(PeriodicSystem& system, Core::Log& log) = 0;

  /**
   * @brief The main functionality of the UnitCellGeometryOptimizer.
   *
   * This function wraps the optimize functions of the underlying optimizers.
   * It has several overloads which take the periodic boundaries from different
   * data structures.
   *
   * @param atoms The AtomCollection to be optimized.
   * @param pbc The PeriodicBoundaries to be optimized.
   * @param log The logger to which eventual output is written.
   * @return int  The final number of optimization cycles carried out.
   */
  virtual int optimize(AtomCollection& atoms, PeriodicBoundaries& pbc, Core::Log& log) = 0;

  // @brief which cell length shall be optimized
  std::vector<bool> optLengths = {true, true, true};
  // @brief if cell angles shall be optimized
  bool optAngles = true;
  // @brief the maximum number of geometry micro iterations
  int geoMaxIterations = 100;
  // @brief the maximum number of unitcell micro iterations
  int cellMaxIterations = 100;
};
/**
 * @brief Settings for a UnitCellGeometryOptimizer
 *
 * Uses template arguments in order to automatically include the
 * settings of underlying objects into the given settings.
 *
 * @tparam OptimizerTypeGeometry The underlying Optimizer class for the geometry.
 * @tparam OptimizerTypeCell The underlying Optimizer class for the unitcell.
 * @tparam ConvergenceCheckType The underlying ConvergenceCheck class.
 */
template<class OptimizerTypeGeometry, class OptimizerTypeCell, class ConvergenceCheckType>
class UnitCellGeometryOptimizerSettings : public Settings {
 public:
  /**
   * @brief Construct a new UnitCellGeometryOptimizerSettings object.
   *
   * Sets the default values of the settings to the current values set in the objects
   * given to the constructor.
   *
   * @param base The base optimizer class.
   * @param geoOptimizer The optimizer for the geometry.
   * @param cellOptimizer The optimizer for the unitcell.
   * @param check The convergence check criteria.
   */
  UnitCellGeometryOptimizerSettings(const UnitCellGeometryOptimizerBase& base,
                                    const GeometryOptimizer<OptimizerTypeGeometry> geoOptimizer,
                                    const OptimizerTypeCell& cellOptimizer, const ConvergenceCheckType& check)
    : Settings("UnitCellGeometryOptimizerSettings") {
    // populate with fields of GeometryOptimizer
    auto geometryOptimizerSettings =
        GeometryOptimizerSettings<OptimizerTypeGeometry, ConvergenceCheckType>(base, geoOptimizer.optimizer, check);
    auto fields = geometryOptimizerSettings.getDescriptorCollection();
    for (const auto& field : fields) {
      this->_fields.push_back(field.first, field.second);
    }
    // add those of cell optimizer if it is different
    if (!std::is_same<OptimizerTypeCell, OptimizerTypeGeometry>::value) {
      cellOptimizer.addSettingsDescriptors(this->_fields);
    }

    // custom convergence for micro iterations
    UniversalSettings::IntDescriptor geoopt_max_iterations(
        "Set the max iterations for the micro geometry optimization cycles.");
    geoopt_max_iterations.setDefaultValue(base.geoMaxIterations);
    this->_fields.push_back(UnitCellGeometryOptimizerBase::geooptMaxIterations, std::move(geoopt_max_iterations));

    UniversalSettings::IntDescriptor cellopt_max_iterations(
        "Set the max iterations for the micro cell optimization cycles.");
    cellopt_max_iterations.setDefaultValue(base.cellMaxIterations);
    this->_fields.push_back(UnitCellGeometryOptimizerBase::celloptMaxIterations, std::move(cellopt_max_iterations));

    // angles
    UniversalSettings::BoolDescriptor angle_opt("Set true if you want to optimize the angles of the unitcell.");
    angle_opt.setDefaultValue(base.optAngles);
    this->_fields.push_back(UnitCellGeometryOptimizerBase::optimizeAnglesKey, std::move(angle_opt));

    // lengths
    UniversalSettings::BoolDescriptor a_opt("Set true if you want to optimize the length of vector a.");
    a_opt.setDefaultValue(base.optLengths[0]);
    UniversalSettings::BoolDescriptor b_opt("Set true if you want to optimize the length of vector b.");
    b_opt.setDefaultValue(base.optLengths[1]);
    UniversalSettings::BoolDescriptor c_opt("Set true if you want to optimize the length of vector c.");
    c_opt.setDefaultValue(base.optLengths[2]);
    this->_fields.push_back(UnitCellGeometryOptimizerBase::optimizeAKey, std::move(a_opt));
    this->_fields.push_back(UnitCellGeometryOptimizerBase::optimizeBKey, std::move(b_opt));
    this->_fields.push_back(UnitCellGeometryOptimizerBase::optimizeCKey, std::move(c_opt));

    this->resetToDefaults();
  }
};

/**
 * @brief Basically just a templated version of the base class UnitCellGeometryOptimizerBase,
 *        where the template defines the actual optimizer used in the optimization.
 *
 * @tparam OptimizerTypeGeometry Expects any of the Optimizer classes. Note that some special optimizers
 *                               may not yet be supported or may need additional specialization.
 * @tparam OptimizerTypeCell Expects any of the Optimizer classes. Note that some special optimizers
 *                           may not yet be supported or may need additional specialization.
 */
template<class OptimizerTypeGeometry, class OptimizerTypeCell>
class UnitCellGeometryOptimizer : public UnitCellGeometryOptimizerBase {
 public:
  /**
   * @brief Construct a new UnitCellGeometryOptimizer object.
   * @param calculator The calculator to be used for the single point/gradient/stresstensor calculations.
   */
  explicit UnitCellGeometryOptimizer(Core::Calculator& calculator) : geoOptimizer(calculator), _calculator(calculator) {
    // WARNING: check for null reference if accessing calculator here, because ReaDuct can plug in a nullptr.
    if (std::is_same<OptimizerTypeCell, NewtonRaphson>::value || std::is_same<OptimizerTypeCell, Bofill>::value ||
        std::is_same<OptimizerTypeCell, EigenvectorFollowing>::value) {
      throw std::logic_error("UnitCellGeometryOptimizer currently does not support Hessian based calculators for the "
                             "cell optimization part.");
    }
    if (std::is_same<OptimizerTypeCell, Dimer>::value || std::is_same<OptimizerTypeCell, Bofill>::value ||
        std::is_same<OptimizerTypeCell, EigenvectorFollowing>::value) {
      throw std::logic_error("UnitCellGeometryOptimizer currently does not support transition state optimizers for "
                             "cell optimization part.");
    }
  };

  /**
   * @brief See UnitCellGeometryOptimizerBase::optimize().
   * @param atoms The AtomCollection to be optimized.
   * @note The unit cell is taken from the calculator.
   * @param log The logger to which eventual output is written.
   *
   * @return int  The final number of optimization cycles carried out.
   */
  int optimize(AtomCollection& atoms, Core::Log& log) final {
    // we are constructing our periodic system from the calculator
    if (!_calculator.settings().valueExists(SettingsNames::periodicBoundaries)) {
      throw std::logic_error("UnitCellGeometryOptimizer requires a calculator that has '" +
                             std::string(SettingsNames::periodicBoundaries) + "' as a possible setting.");
    }
    auto pbc = Utils::PeriodicBoundaries(_calculator.settings().getString(SettingsNames::periodicBoundaries));
    int cycles = optimize(atoms, pbc, log);
    // be sure that the calculator is up-to-date
    _calculator.settings().modifyString(SettingsNames::periodicBoundaries, pbc.getPeriodicBoundariesString());
    return cycles;
  }

  /**
   * @brief See UnitCellGeometryOptimizerBase::optimize().
   * @param atoms The AtomCollection to be optimized.
   * @param pbc The unit cell to be optimized.
   * @param log The logger to which eventual output is written.
   *
   * @return int  The final number of optimization cycles carried out.
   */
  int optimize(AtomCollection& atoms, PeriodicBoundaries& pbc, Core::Log& log) final {
    auto system = Utils::PeriodicSystem(pbc, atoms);
    int cycles = optimize(system, log);
    atoms = system.atoms;
    pbc = system.pbc;
    return cycles;
  }
  /**
   * @brief See UnitCellGeometryOptimizerBase::optimize().
   * @param system The Periodicsystem to be optimized.
   * @param log The logger to which eventual output is written.
   *
   * @return int  The final number of optimization cycles carried out.
   */
  int optimize(PeriodicSystem& system, Core::Log& log) final {
    if (!_calculator.possibleProperties().containsSubSet(Utils::Property::StressTensor)) {
      throw std::logic_error("UnitCellGeometryOptimizer requires a calculator that can calculate a stress tensor.");
    }
    if (!_calculator.settings().valueExists(SettingsNames::periodicBoundaries)) {
      throw std::logic_error("UnitCellGeometryOptimizer requires a calculator that has '" +
                             std::string(SettingsNames::periodicBoundaries) + "' as a possible setting.");
    }
    // ensure our settings are up to date in case members have been altered + sanity checks
    setSettings(getSettings());
    bool disabledCellOpt = !optAngles && std::find(optLengths.begin(), optLengths.end(), true) == optLengths.end();
    if (disabledCellOpt) {
      log.warning
          << "Angle and all length optimizations have been disabled. The following optimization solely optimizes "
             "the structure and not the unitcell."
          << Core::Log::endl;
    }
    _calculator.setStructure(system.atoms);
    _calculator.settings().modifyString(SettingsNames::periodicBoundaries, system.pbc.getPeriodicBoundariesString());
    _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::StressTensor);
    unsigned int totalCycles = 0;
    bool converged = false;
    while (!converged) {
      int cellCycles = 0;
      if (!disabledCellOpt) {
        cellCycles = optimizeCell(system, log);
        log.output << "    Number of unitcell micro optimization steps: " << std::to_string(cellCycles) << Core::Log::nl;
        totalCycles += cellCycles;
        if (totalCycles >= check.maxIter) {
          log.warning << "Reached total max iterations during unitcell geometry optimization." << Core::Log::endl;
          return totalCycles;
        }
      }
      int geoCycles = geoOptimizer.optimize(system.atoms, log);
      log.output << "    Number of geometry micro optimization steps: " << std::to_string(geoCycles) << Core::Log::nl;
      totalCycles += geoCycles;
      if (totalCycles >= check.maxIter) {
        log.warning << "Reached total max iterations during unitcell geometry optimization." << Core::Log::endl;
        return totalCycles;
      }
      // BFGS per default does not immediately check for convergence -> change this after first geometry/cell cycle
      // have to go long way with settings instead of directly altering members, because compilation would fail for
      // optimizers that are missing the 'minIter' member because of template deduction
      if (!disabledCellOpt && instanceof <Bfgs>(geoOptimizer.optimizer)) {
        auto settings = this->getSettings();
        settings.modifyInt(Bfgs::bfgsMinIter, 1);
        this->setSettings(settings);
      }
      if (!disabledCellOpt && instanceof <Bfgs>(cellOptimizer)) {
        auto settings = this->getSettings();
        settings.modifyInt(Bfgs::bfgsMinIter, 1);
        this->setSettings(settings);
      }
      // we are converged if both optimizer did not need to change their given parameters
      // or we do not optimize cell and geometry optimization was below maximum cycles
      converged = (cellCycles == 2 && geoCycles == 2) || (disabledCellOpt && geoCycles < geoMaxIterations);
    }
    return totalCycles;
  }

  /**
   * @brief Function to apply the given settings to underlying classes.
   * @param settings The new settings.
   */
  void setSettings(const Settings& settings) {
    geoOptimizer.setSettings(settings);
    check.applySettings(settings);
    cellOptimizer.applySettings(settings);

    // For Cartesian constraints:
    this->fixedAtoms = settings.getIntList(GeometryOptimizerBase::geooptFixedAtomsKey);
    if (!this->fixedAtoms.empty()) {
      throw std::logic_error("Cartesian constraints cannot be set when the unit cell shall be optimized");
    }

    // convergence settings
    this->geoMaxIterations = settings.getInt(UnitCellGeometryOptimizerBase::geooptMaxIterations);
    Settings s1 = geoOptimizer.getSettings();
    s1.modifyInt(GradientBasedCheck::gconvMaxIterKey, geoMaxIterations);
    geoOptimizer.setSettings(s1);
    // celloptimer is plain optimizer and we work around this in private function before optimize call
    this->cellMaxIterations = settings.getInt(UnitCellGeometryOptimizerBase::celloptMaxIterations);

    // unit cell optimization settings
    this->optAngles = settings.getBool(UnitCellGeometryOptimizerBase::optimizeAnglesKey);
    bool a = settings.getBool(UnitCellGeometryOptimizerBase::optimizeAKey);
    bool b = settings.getBool(UnitCellGeometryOptimizerBase::optimizeBKey);
    bool c = settings.getBool(UnitCellGeometryOptimizerBase::optimizeCKey);
    this->optLengths = {a, b, c};
  };
  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  Settings getSettings() const {
    return UnitCellGeometryOptimizerSettings<OptimizerTypeGeometry, OptimizerTypeCell, GradientBasedCheck>(
        *this, geoOptimizer, cellOptimizer, check);
  };
  /**
   * @brief Get a copy of the settings of the calculator used for the energy calculations during the optimization.
   * @return std::shared_ptr<Settings> The settings of the calculator.
   */
  std::shared_ptr<Settings> getCalculatorSettings() const {
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
  void addObserver(std::function<void(const int&, const double&, const Eigen::VectorXd&)> function) {
    geoOptimizer.addObserver(function);
    cellOptimizer.addObserver(function);
  }
  /**
   * @brief Clear all existing observer functions.
   *
   * For optimization problems with very fast evaluations of the underlying function
   * the removal of all observers can increase performance as the observers are given
   * as std::functions and can not be added via templates.
   */
  void clearObservers() {
    geoOptimizer.clearObservers();
    cellOptimizer.clearObservers();
  }
  /**
   * @brief The underlying convergence check
   *
   * @note getter to be accessible via base class
   * @return GradientBasedCheck the class holding all convergence thresholds.
   */
  const GradientBasedCheck& getConvergenceCheck() const {
    return check;
  };
  ~UnitCellGeometryOptimizer() override = default;
  /// @brief The underlying GeometryOptimizer, public in order to change its settings.
  GeometryOptimizer<OptimizerTypeGeometry> geoOptimizer;
  /// @brief The underlying optimizer, public in order to change its settings.
  OptimizerTypeCell cellOptimizer;
  /// @brief The underlying convergence check, public in order to change its settings.
  GradientBasedCheck check;

 private:
  Core::Calculator& _calculator;

  /**
   * @brief Utility to flatten derivative matrix to vector and apply constraints.
   *
   * @param latticeDerivative The derivative in matrix form without constraints.
   * @param pbc The unit cell.
   * @return Eigen::VectorXd The flattened and constraint gradient vector for unit cell opt step.
   */
  inline Eigen::VectorXd flattenLatticeDerivativeAndApplyConstraints(Eigen::Matrix3d latticeDerivative,
                                                                     const PeriodicBoundaries& pbc) {
    Eigen::VectorXd result = Eigen::Map<Eigen::VectorXd>(latticeDerivative.data(), latticeDerivative.size());
    if (!optAngles) { // angles are constrained
                      // project forces to be along original lattice
      latticeDerivative = latticeDerivative * pbc.getCellMatrixWithNormalizedVectors();
      result = Eigen::VectorXd::Zero(9);
      for (int i = 0; i < 3; ++i) {
        if (optLengths[i]) {
          result[3 * i + i] += latticeDerivative.col(i).sum();
        }
      }
    }
    // nullify constrained lengths
    for (int i = 0; i < 3; ++i) {
      if (!optLengths[i]) {
        result[3 * i + i] = 0.0;
      }
    }
    return result;
  }

  // @brief utility overload
  inline int optimizeCell(PeriodicSystem& system, Core::Log& log) {
    return optimizeCell(system.atoms, system.pbc, log);
  }

  /**
   * @brief Equivalent of optimize of GeometryOptimizer but solely for unit cell optimization.
   *
   * @param atoms The atoms in the unit cell, which are slightly shifted with the unit cell.
   * @param pbc The unit cell.
   * @param log The logger for the optimization.
   * @return Eigen::VectorXd The flattened and constraint gradient vector for unit cell opt step.
   */
  inline int optimizeCell(AtomCollection& atoms, Utils::PeriodicBoundaries& pbc, Core::Log& log) {
    // we first define our update function taking parameters as vector (flattened unit cell definition)
    // and updating energy and gradients
    auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
      // set new PBCs from parameters and modify atom positions via relative coordinates
      Eigen::Matrix3d cellMatrix = Eigen::Map<const Eigen::Matrix3d>(parameters.data(), 3, 3);
      PositionCollection relCoordinates = pbc.transform(atoms.getPositions(), false);
      pbc.setCellMatrix(cellMatrix);
      PositionCollection coordinates = pbc.transform(relCoordinates);
      atoms.setPositions(coordinates);
      // new calculation
      _calculator.modifyPositions(coordinates);
      _calculator.setRequiredProperties(Property::Energy | Property::Gradients | Property::StressTensor);
      _calculator.settings().modifyString(SettingsNames::periodicBoundaries, pbc.getPeriodicBoundariesString());
      Results results =
          CalculationRoutines::calculateWithCatch(_calculator, log, "Aborting optimization due to failed calculation");
      value = results.get<Property::Energy>();
      // set gradient
      Eigen::Matrix3d stressTensor = results.get<Property::StressTensor>();
      Eigen::Matrix3d latticeDerivative = stressTensor * pbc.getInverseCellMatrix();
      latticeDerivative *= (-pbc.getVolume()); // using DFTB+ sign convention here
      gradients = flattenLatticeDerivativeAndApplyConstraints(latticeDerivative, pbc);
    };
    // transform cell matrix, which are the parameters to be optimized, to vector
    Eigen::VectorXd positions = Eigen::Map<Eigen::VectorXd>(pbc.getCellMatrix().data(), pbc.getCellMatrix().size());
    // hack for keeping maximum cycles
    auto maxTotalCycles = check.maxIter;
    check.maxIter = static_cast<unsigned>(this->cellMaxIterations);
    // optimize
    int cycles = cellOptimizer.optimize(positions, update, check, log);
    check.maxIter = maxTotalCycles;
    return cycles;
  }

  // @brief utility to check if we got a specific optimizer
  // reason is that std::is_same still tries out code in compilation
  template<typename Base, typename T>
  inline bool instanceof (const T) {
    return std::is_base_of<Base, T>::value;
  }
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_UNITCELLGEOMETRYOPTIMIZER_H_
