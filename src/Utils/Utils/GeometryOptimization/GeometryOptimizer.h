/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_GEOMETRYOPTIMIZER_H_
#define UTILS_GEOMETRYOPTIMIZER_H_

#include "Utils/CalculatorBasics/PropertyList.h"
#include "Utils/CalculatorBasics/Results.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/InternalCoordinates.h"
#include "Utils/Optimizer/GradientBased/Bfgs.h"
#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/GradientBased/Lbfgs.h"
#include "Utils/Optimizer/HessianBased/Bofill.h"
#include "Utils/Optimizer/HessianBased/EigenvectorFollowing.h"
#include "Utils/Optimizer/HessianBased/NewtonRaphson.h"
#include <Core/Interfaces/Calculator.h>
#include <Eigen/Core>
#include <iostream>

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
  static constexpr const char* geooptTransfromCoordinatesKey = "geoopt_transform_coordinates";
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
   * @return int  The final number of optimization cycles carried out.
   */
  virtual int optimize(AtomCollection& atoms) = 0;
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
   * @brief Switch to transform the coordinates from Cartesian into an internal space.
   *
   * The optimization will be carried out in the internal coordinate space possibly
   * accellerating convergence.
   */
  bool transformCoordinates = true;
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

    UniversalSettings::BoolDescriptor geoopt_transform_coordinates(
        "Switch to transform the coordinates from Cartesian into an internal space.");
    geoopt_transform_coordinates.setDefaultValue(base.transformCoordinates);
    this->_fields.push_back(GeometryOptimizerBase::geooptTransfromCoordinatesKey, geoopt_transform_coordinates);
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
  GeometryOptimizer(Core::Calculator& calculator) : _calculator(calculator){};
  /**
   * @brief See GeometryOptimizerBase::optimize().
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
    const unsigned int nAtoms = atoms.size();
    auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
      Utils::PositionCollection coordinates;
      if (this->transformCoordinates) {
        coordinates = transformation->coordinatesToCartesian(parameters);
      }
      else {
        coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
      }
      _calculator.modifyPositions(coordinates);
      _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
      atoms.setPositions(coordinates);
      Utils::Results results = _calculator.calculate("Geometry Optimization Cycle");
      value = results.get<Property::Energy>();
      if (this->transformCoordinates) {
        gradients = transformation->gradientsToInternal(results.get<Property::Gradients>());
      }
      else {
        gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
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
    this->transformCoordinates = settings.getBool(GeometryOptimizerBase::geooptTransfromCoordinatesKey);
  };
  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  virtual Settings getSettings() const override {
    return GeometryOptimizerSettings<OptimizerType, GradientBasedCheck>(*this, optimizer, check);
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

/*=============================*
 *  Approximate Hessian Based
 *=============================*/

template<>
inline int GeometryOptimizer<Bfgs>::optimize(AtomCollection& atoms) {
  // Configure Calculator
  _calculator.setStructure(atoms);
  _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  // Transform into internal coordinates
  auto transformation = std::make_shared<InternalCoordinates>(atoms);
  // Define update function
  const unsigned int nAtoms = atoms.size();
  auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
    Utils::PositionCollection coordinates;
    if (this->transformCoordinates) {
      coordinates = transformation->coordinatesToCartesian(parameters);
    }
    else {
      coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
    }
    _calculator.modifyPositions(coordinates);
    _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
    atoms.setPositions(coordinates);
    Utils::Results results = _calculator.calculate("Geometry Optimization Cycle");
    value = results.get<Property::Energy>();
    if (this->transformCoordinates) {
      gradients = transformation->gradientsToInternal(results.get<Property::Gradients>());
    }
    else {
      gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
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
    optimizer.invH = transformation->inverseHessianGuess();
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
inline int GeometryOptimizer<Bofill>::optimize(AtomCollection& atoms) {
  // Configure Calculator
  _calculator.setStructure(atoms);
  _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
  // Transformation into internal basis
  Eigen::MatrixXd transformation;
  auto elements = atoms.getElements();
  if (this->transformCoordinates) {
    transformation = Geometry::calculateRotTransFreeTransformMatrix(atoms.getPositions(), elements);
  }
  // Define update function
  const unsigned int nAtoms = atoms.size();
  auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                          Eigen::MatrixXd& hessian, bool calcHessian) {
    Utils::PositionCollection coordinates;
    if (this->transformCoordinates) {
      auto tmp = (transformation * parameters).eval();
      coordinates = Eigen::Map<const Utils::PositionCollection>(tmp.data(), nAtoms, 3);
    }
    else {
      coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
    }
    _calculator.modifyPositions(coordinates);
    atoms.setPositions(coordinates);

    if (calcHessian) {
      _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
      Utils::Results results = _calculator.calculate("Geometry Optimization Cycle");
      value = results.get<Property::Energy>();
      if (this->transformCoordinates) {
        auto tmp = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
        gradients = (transformation.transpose() * tmp).eval();
        hessian = transformation.transpose() * results.get<Property::Hessian>() * transformation;
      }
      else {
        gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
        hessian = results.get<Property::Hessian>();
      }
    }
    else {
      _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
      Utils::Results results = _calculator.calculate("Geometry Optimization Cycle");
      value = results.get<Property::Energy>();
      if (this->transformCoordinates) {
        auto tmp = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
        gradients = (transformation.transpose() * tmp).eval();
      }
      else {
        gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
      }
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

template<>
inline int GeometryOptimizer<NewtonRaphson>::optimize(AtomCollection& atoms) {
  // Configure Calculator
  _calculator.setStructure(atoms);
  _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
  // Transformation into internal basis
  Eigen::MatrixXd transformation;
  auto elements = atoms.getElements();
  if (this->transformCoordinates) {
    transformation = Geometry::calculateRotTransFreeTransformMatrix(atoms.getPositions(), elements);
  }
  // Define update function
  const unsigned int nAtoms = atoms.size();
  auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                          Eigen::MatrixXd& hessian, bool /*calcHessian*/) {
    Utils::PositionCollection coordinates;
    if (this->transformCoordinates) {
      auto tmp = (transformation * parameters).eval();
      coordinates = Eigen::Map<const Utils::PositionCollection>(tmp.data(), nAtoms, 3);
    }
    else {
      coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
    }
    _calculator.modifyPositions(coordinates);
    atoms.setPositions(coordinates);
    _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
    Utils::Results results = _calculator.calculate("Geometry Optimization Cycle");
    value = results.get<Property::Energy>();
    if (this->transformCoordinates) {
      auto tmp = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
      gradients = (transformation.transpose() * tmp).eval();
      hessian = transformation.transpose() * results.get<Property::Hessian>() * transformation;
    }
    else {
      gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
      hessian = results.get<Property::Hessian>();
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

template<>
inline int GeometryOptimizer<EigenvectorFollowing>::optimize(AtomCollection& atoms) {
  // Configure Calculator
  _calculator.setStructure(atoms);
  _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
  // Transformation into internal basis
  Eigen::MatrixXd transformation;
  auto elements = atoms.getElements();
  if (this->transformCoordinates) {
    transformation = Geometry::calculateRotTransFreeTransformMatrix(atoms.getPositions(), elements);
  }
  // Define update function
  const unsigned int nAtoms = atoms.size();
  auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                          Eigen::MatrixXd& hessian, bool /*calcHessian*/) {
    Utils::PositionCollection coordinates;
    if (this->transformCoordinates) {
      auto tmp = (transformation * parameters).eval();
      coordinates = Eigen::Map<const Utils::PositionCollection>(tmp.data(), nAtoms, 3);
    }
    else {
      coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
    }
    _calculator.modifyPositions(coordinates);
    atoms.setPositions(coordinates);
    _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
    Utils::Results results = _calculator.calculate("Geometry Optimization Cycle");
    value = results.get<Property::Energy>();
    if (this->transformCoordinates) {
      auto tmp = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
      gradients = (transformation.transpose() * tmp).eval();
      hessian = transformation.transpose() * results.get<Property::Hessian>() * transformation;
    }
    else {
      gradients = Eigen::Map<const Eigen::VectorXd>(results.get<Property::Gradients>().data(), nAtoms * 3);
      hessian = results.get<Property::Hessian>();
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

} // namespace Utils
} // namespace Scine

#endif // UTILS_GEOMETRYOPTIMIZER_H_
