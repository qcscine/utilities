/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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
#include "Utils/UniversalSettings/OptimizationSettingsNames.h"
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
    this->_fields.push_back(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem,
                            std::move(geoopt_coordinate_system));

    UniversalSettings::IntListDescriptor geooptFixedAtoms(
        "List of atoms with Cartesian constraints applied to them during the optimization.");
    geooptFixedAtoms.setDefaultValue(base.fixedAtoms);
    this->_fields.push_back(SettingsNames::Optimizations::GeometryOptimizer::fixedAtoms, std::move(geooptFixedAtoms));

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
  explicit GeometryOptimizer(Core::Calculator& calculator) : _calculator(calculator) {
    // WARNING: check for null reference if accessing calculator here, because ReaDuct can plug in a nullptr.
    /* set private members according to template */
    // Performs better with Cartesian coordinates than with current internal coordinates
    if (std::is_same<OptimizerType, Dimer>::value) {
      this->coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
    }
    // With current internal coordinates not possible
    if (std::is_same<OptimizerType, Lbfgs>::value || std::is_same<OptimizerType, NewtonRaphson>::value ||
        std::is_same<OptimizerType, Bofill>::value || std::is_same<OptimizerType, EigenvectorFollowing>::value) {
      // TODO: fix the Hessian projection in the L-BFGS to allow for internals, others require Hessian in internals
      this->coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
      this->_internalAvailable = false;
    }
  };
  /**
   * @brief See GeometryOptimizerBase::optimize().
   * @param atoms The AtomCollection (Geometry) to be optimized.
   * @param log The logger to which eventual output is written.
   *
   * @return int  The final number of optimization cycles carried out.
   */
  int optimize(AtomCollection& atoms, Core::Log& log) final {
    // Configure members
    _atoms = atoms;
    _log = std::make_shared<Core::Log>(log);
    auto calcStructure = _calculator.getStructure();
    Results oldResults = _calculator.results();
    if (calcStructure && calcStructure->getElements() == atoms.getElements()) {
      const Eigen::MatrixXd& oldPositions = _calculator.getPositions();
      _calculator.modifyPositions(atoms.getPositions());
      // Check in Cartesian Coordinates if the given positions and the positions in the calculator are identical
      if (oldPositions.rows() == atoms.getPositions().rows() && atoms.getPositions().isApprox(oldPositions, 1e-12)) {
        _calculator.results() = oldResults;
        _useOldResults = true;
      }
    }
    else {
      _calculator.setStructure(atoms);
    }
    auto originalProperties = _calculator.getRequiredProperties();
    originalProperties.addProperties(_requiredProperties);
    _calculator.setRequiredProperties(originalProperties);
    // Transformation into internal coordinates
    if (this->coordinateSystem == CoordinateSystem::Internal) {
      if (!_internalAvailable) {
        throw std::logic_error("Internal coordinates are currently not supported for this optimizer.");
      }
      _transformation = std::make_shared<InternalCoordinates>(atoms);
    }
    else if (this->coordinateSystem == CoordinateSystem::CartesianWithoutRotTrans) {
      _transformation = std::make_shared<InternalCoordinates>(atoms, true);
    }
    // Get initial positions
    const unsigned int nAtoms = atoms.size();
    Eigen::VectorXd positions;
    if (_transformation) {
      positions = _transformation->coordinatesToInternal(atoms.getPositions());
    }
    else {
      positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), nAtoms * 3);
      // set constraints
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
      cycles = callOptimizer(positions);
    }
    catch (const InternalCoordinatesException& e) {
      log.output << "Internal coordinates broke down. Continuing in Cartesians." << Core::Log::nl;
      // Update coordinates to the last ones that were successfully reconverted
      Utils::PositionCollection lastCoordinates = _calculator.getPositions();
      // Disable true internal coordinates
      this->coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
      _transformation = std::make_shared<InternalCoordinates>(_atoms.get(), true);
      Eigen::VectorXd lastPositions = _transformation->coordinatesToInternal(lastCoordinates);
      // Restart optimization
      optimizer.prepareRestart(optimizer.getCycle());
      cycles = callOptimizer(lastPositions);
      positions = lastPositions;
    }
    // Update Atom collection and return
    Utils::PositionCollection coordinates;
    if (_transformation) {
      coordinates = _transformation->coordinatesToCartesian(positions);
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
        settings.getString(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem));

    // For Cartesian constraints:
    this->fixedAtoms = settings.getIntList(SettingsNames::Optimizations::GeometryOptimizer::fixedAtoms);

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
  /// @brief The underlying optimizer, public in order to change its settings.
  OptimizerType optimizer;
  /// @brief The underlying convergence check, public in order to change its settings.
  GradientBasedCheck check;

 private:
  Core::Calculator& _calculator;
  std::shared_ptr<InternalCoordinates> _transformation;
  std::shared_ptr<Core::Log> _log;
  // we need to store the atoms as class members to manipulate them in a general lambda function outside of the given
  // atom scope; they need to be stored such that their manipulation also manipulates the given atoms
  // used reference_wrapper, which gets a dummy init value to make the class constructable without atoms
  AtomCollection _dummy = AtomCollection();
  std::reference_wrapper<AtomCollection> _atoms = std::ref(_dummy);
  // following members may be adapted in constructur depending on template
  PropertyList _requiredProperties = Utils::Property::Energy | Utils::Property::Gradients;
  bool _internalAvailable = true;
  bool _useOldResults = false;

  /**
   * @brief Lambda function passed to mathematical optimizer, which calls it to update value and gradient
   */
  const std::function<void(const Eigen::VectorXd&, double&, Eigen::VectorXd&)> _gradientUpdateLambda =
      [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
        const unsigned int nAtoms = _atoms.get().size();
        Utils::PositionCollection coordinates;
        if (_transformation) {
          coordinates = _transformation->coordinatesToCartesian(parameters);
        }
        else {
          coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
        }
        _calculator.modifyPositions(coordinates);
        auto originalProperties = _calculator.getRequiredProperties();
        if (originalProperties.containsSubSet(Property::Hessian)) {
          originalProperties.removeProperty(Property::Hessian);
        }
        if (originalProperties.containsSubSet(Property::PartialHessian)) {
          originalProperties.removeProperty(Property::PartialHessian);
        }
        if (originalProperties.containsSubSet(Property::Thermochemistry)) {
          originalProperties.removeProperty(Property::Thermochemistry);
        }
        originalProperties.addProperties(_requiredProperties);
        _calculator.setRequiredProperties(originalProperties);
        _atoms.get().setPositions(coordinates);
        Results results =
            CalculationRoutines::calculateWithCatch(_calculator, *_log, "Aborting optimization due to failed calculation");
        value = results.get<Property::Energy>();
        if (_transformation) {
          gradients = _transformation->gradientsToInternal(results.get<Property::Gradients>());
        }
        else {
          gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
        }
      };
  /**
   * @brief Lambda function passed to mathematical optimizer, which calls it to update value, gradient, and optionally
   * Hessian
   */
  const std::function<void(const Eigen::VectorXd&, double&, Eigen::VectorXd&, Eigen::MatrixXd&, bool)> _hessianUpdateLambda =
      [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients, Eigen::MatrixXd& hessian,
          bool calcHessian) {
        if (!calcHessian) {
          _gradientUpdateLambda(parameters, value, gradients);
        }
        else {
          const unsigned int nAtoms = _atoms.get().size();
          Utils::PositionCollection coordinates;
          if (_transformation) {
            coordinates = _transformation->coordinatesToCartesian(parameters);
          }
          else {
            coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
          }
          Results results = _calculator.results();
          bool recalculate = true;

          if (_useOldResults) {
            _useOldResults = false;
            if (results.has<Utils::Property::Energy>() && results.has<Utils::Property::Gradients>() &&
                results.has<Utils::Property::Hessian>()) {
              // Same coordinates and all results are present, no need to calculate anything
              recalculate = false;
            }
          }
          if (recalculate) {
            _calculator.modifyPositions(coordinates);
            auto originalProperties = _calculator.getRequiredProperties();
            originalProperties.addProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
            _calculator.setRequiredProperties(originalProperties);
            results = CalculationRoutines::calculateWithCatch(_calculator, *_log,
                                                              "Aborting optimization due to failed calculation");
          }
          _atoms.get().setPositions(coordinates);
          value = results.get<Property::Energy>();
          if (_transformation) {
            gradients = _transformation->gradientsToInternal(results.get<Property::Gradients>());
            hessian = _transformation->hessianToInternal(results.get<Property::Hessian>());
          }
          else {
            gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
            hessian = results.get<Property::Hessian>();
          }
        }
      };

  /**
   * @brief Prepare underlying optimizer, default is to do nothing, since most don't need preparation
   */
  inline void optimizerPreparation() {
    return;
  };
  /**
   * @brief Call the underlying optimizer with a lambda update function; default is gradient update function
   * @param positions The coordinates of the atoms as a vector, which are optimized
   * @return int the number of required cycles
   */
  inline int callOptimizer(Eigen::VectorXd& positions) {
    optimizerPreparation();
    return optimizer.optimize(positions, _gradientUpdateLambda, check, *_log);
  };
};

/*=====================================================================================*
 *  Preparation specification via template, set invH and transformation for internals
 *=====================================================================================*/
template<>
inline void GeometryOptimizer<Bfgs>::optimizerPreparation() {
  if (this->coordinateSystem == CoordinateSystem::Internal) {
    optimizer.invH = _transformation->inverseHessianGuess();
    // pass a weak ptr to lambda function as shared ptr may cause leaks:
    // https://floating.io/2017/07/lambda-shared_ptr-memory-leak/
    std::weak_ptr<InternalCoordinates> transformation(_transformation);
    optimizer.projection = [transformation](Eigen::MatrixXd& inv) {
      auto t = transformation.lock();
      inv = t->projectHessianInverse(inv);
    };
  }
  else {
    optimizer.invH.resize(0, 0);
    optimizer.projection = nullptr;
  }
}
template<>
inline void GeometryOptimizer<Dimer>::optimizerPreparation() {
  if (this->coordinateSystem == CoordinateSystem::Internal) {
    optimizer.invH = _transformation->inverseHessianGuess();
    // pass a weak ptr to lambda function as shared ptr may cause leaks:
    // https://floating.io/2017/07/lambda-shared_ptr-memory-leak/
    std::weak_ptr<InternalCoordinates> transformation(_transformation);
    optimizer.projection = std::make_shared<std::function<void(Eigen::MatrixXd&)>>([transformation](Eigen::MatrixXd& inv) {
      auto t = transformation.lock();
      inv = t->projectHessianInverse(inv);
    });
  }
  else {
    optimizer.projection = nullptr;
  }
}

/*====================================================================================*
 *  Update specification via template, overwrite gradient update with Hessian update
 *====================================================================================*/
template<>
inline int GeometryOptimizer<Bofill>::callOptimizer(Eigen::VectorXd& positions) {
  optimizerPreparation();
  return optimizer.optimize(positions, _hessianUpdateLambda, check, *_log);
}
template<>
inline int GeometryOptimizer<EigenvectorFollowing>::callOptimizer(Eigen::VectorXd& positions) {
  optimizerPreparation();
  return optimizer.optimize(positions, _hessianUpdateLambda, check, *_log);
}
template<>
inline int GeometryOptimizer<NewtonRaphson>::callOptimizer(Eigen::VectorXd& positions) {
  optimizerPreparation();
  return optimizer.optimize(positions, _hessianUpdateLambda, check, *_log);
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_GEOMETRYOPTIMIZER_H_
