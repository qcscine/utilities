/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/GeometryOptimization/NtOptimizer.h"
#include "Utils/CalculatorBasics/CalculationRoutines.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/GeometryOptimization/NtOptimizerSettings.h"
#include "Utils/Optimizer/GradientBased/Bfgs.h"
#include <array>
#include <boost/exception/diagnostic_information.hpp>
#include <cmath>
#include <valarray>

namespace Scine {
namespace Utils {

int NtOptimizer::optimize(AtomCollection& atoms, Core::Log& log) {
  this->sanityCheck(atoms);
  // Configure Calculator
  _calculator.setStructure(atoms);
  _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  Utils::PositionCollection coordinates = atoms.getPositions();
  const unsigned int nAtoms = atoms.size();
  int cycle = 0;
  _values.clear();
  for (unsigned int loop = 0; loop < this->check.maxIter; loop++) {
    cycle++;
    // Micro cycles performed in true Cartesians due to constraints
    if (cycle > 1 && useMicroCycles && nAtoms > 2) {
      // Define micro iteration optimizer
      Bfgs bfgs;
      bfgs.projection = nullptr;
      bfgs.trustRadius = 0.2;
      bfgs.useTrustRadius = true;
      // Define micro iteration convergence
      GradientBasedCheck microIterCheck;
      microIterCheck.maxIter = (fixedNumberOfMicroCycles ? numberOfMicroCycles : std::min(cycle, numberOfMicroCycles));
      // Define micro iteration update
      auto const microIterUpdate = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
        coordinates = Eigen::Map<const Utils::PositionCollection>(parameters.data(), nAtoms, 3);
        _calculator.modifyPositions(coordinates);
        _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
        atoms.setPositions(coordinates);
        Results results = CalculationRoutines::calculateWithCatch(_calculator, log, "Calculation in NT optimization failed.");
        value = results.get<Property::Energy>();
        // Apply Cartesian constraints
        auto gradientMatrix = results.get<Property::Gradients>();
        this->updateGradients(atoms, value, gradientMatrix);
        gradients = Eigen::Map<const Eigen::VectorXd>(gradientMatrix.data(), nAtoms * 3);
      };
      // Run micro iterations
      Eigen::VectorXd positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), atoms.size() * 3);
      try {
        bfgs.optimize(positions, microIterUpdate, microIterCheck);
      }
      catch (...) {
        if (cycle > 5) {
          coordinates = this->extractTsGuess();
          atoms.setPositions(coordinates);
          _calculator.modifyPositions(coordinates);
          return cycle;
        }
        throw std::runtime_error("NT micro iterations failed.");
      }
      coordinates = Eigen::Map<const Utils::PositionCollection>(positions.data(), nAtoms, 3);
      atoms.setPositions(coordinates);
    }

    _calculator.modifyPositions(coordinates);

    // Calculate data
    Utils::Results results;
    try {
      results = _calculator.calculate("NT Macro Cycle"); // might throw exception
      if (!results.get<Property::SuccessfulCalculation>()) {
        throw std::runtime_error("Calculator signalled unsuccessful calculation.");
      }
    }
    catch (...) {
      // try extraction if we already walked a bit, calculation failure might occur because atoms are already quite close
      if (cycle > 5) {
        coordinates = this->extractTsGuess();
        atoms.setPositions(coordinates);
        _calculator.modifyPositions(coordinates);
        return cycle;
      }
      throw Core::UnsuccessfulCalculationException("Calculation in NT optimization failed: " +
                                                   boost::current_exception_diagnostic_information());
    }
    double value = results.get<Property::Energy>();
    auto gradients = results.get<Property::Gradients>();
    this->triggerObservers(cycle, value, Eigen::Map<const Eigen::VectorXd>(coordinates.data(), nAtoms * 3));
    // Evaluate additional force
    this->updateGradients(atoms, value, gradients, true);
    _values.push_back(value);
    _trajectory.push_back(coordinates);
    // Check convergence
    if (this->convergedOptimization(atoms)) {
      coordinates = this->extractTsGuess();
      atoms.setPositions(coordinates);
      _calculator.modifyPositions(coordinates);
      return cycle;
    }
    // Update positions (SD)
    this->updateCoordinates(coordinates, atoms, gradients);
    atoms.setPositions(coordinates);
    _calculator.modifyPositions(coordinates);
  }
  return cycle;
}

void NtOptimizer::sanityCheck(const AtomCollection& atoms) const {
  if (lhsList.empty() || rhsList.empty()) {
    throw std::logic_error("Called NT optimization without atoms to apply a force to.");
  }
  // check for negative indices
  if (std::any_of(lhsList.begin(), lhsList.end(), [&](int i) { return i < 0; })) {
    throw std::logic_error(
        "At least one of the given indices in the 'nt_lhs_list' is negative, which is not possible.");
  }
  if (std::any_of(rhsList.begin(), rhsList.end(), [&](int i) { return i < 0; })) {
    throw std::logic_error(
        "At least one of the given indices in the 'nt_rhs_list' is negative, which is not possible.");
  }
  // check for too large indices
  const int nAtoms = atoms.size();
  for (const auto i : lhsList) {
    if (i >= nAtoms) {
      throw std::logic_error("Index '" + std::to_string(i) + "' in 'nt_lhs_list' too large for " +
                             std::to_string(nAtoms) + " atoms.");
    }
  }
  for (const auto i : rhsList) {
    if (i >= nAtoms) {
      throw std::logic_error("Index '" + std::to_string(i) + "' in 'nt_rhs_list' too large for " +
                             std::to_string(nAtoms) + " atoms.");
    }
  }
  // check if there is an identical index on both sides
  for (const auto lhs : lhsList) {
    if (std::find(rhsList.begin(), rhsList.end(), lhs) != rhsList.end()) {
      throw std::logic_error("The index " + std::to_string(lhs) + " is present in lhs and rhs, which is not valid.");
    }
  }
  // check if there is an overlap of reaction indices and constrained atoms and too large indices of restrained atoms
  if (!this->fixedAtoms.empty()) {
    for (const auto& a : fixedAtoms) {
      if (a >= nAtoms) {
        throw std::logic_error("Constrained atom index " + std::to_string(a) + " is too large for " +
                               std::to_string(nAtoms) + " atoms.");
      }
      std::vector<int> movableReactionIndicesList;
      if (this->movableSide == "both") {
        movableReactionIndicesList = lhsList;
        movableReactionIndicesList.insert(movableReactionIndicesList.end(), rhsList.begin(), rhsList.end());
      }
      else if (this->movableSide == "lhs") {
        movableReactionIndicesList = lhsList;
      }
      else if (this->movableSide == "rhs") {
        movableReactionIndicesList = rhsList;
      }
      else {
        // should be handled by settings
        throw std::logic_error("Unknown input " + movableSide + " for " + NtOptimizer::ntMovableSide);
      }
      if (std::find(movableReactionIndicesList.begin(), movableReactionIndicesList.end(), a) !=
          movableReactionIndicesList.end()) {
        throw std::logic_error("The atom index " + std::to_string(a) +
                               " was specified to be constrained although "
                               "it is within a movable reaction side. You can change '" +
                               NtOptimizer::ntMovableSide + "' if you want to constrain atoms.");
      }
    }
  }
}

void NtOptimizer::updateGradients(const AtomCollection& atoms, const double& /* energy */,
                                  GradientCollection& gradients, bool addForce) const {
  // define reaction coordinate based on LHS and RHS list
  GradientCollection reactionCoordinate(gradients);
  reactionCoordinate.setZero();
  Displacement c2c = this->centerToCenterVector(atoms.getPositions());
  c2c /= c2c.norm();
  for (const auto l : lhsList) {
    reactionCoordinate.row(l) = -c2c;
  }
  for (const auto r : rhsList) {
    reactionCoordinate.row(r) = c2c;
  }
  /* Update gradients */
  // remove gradient along reaction coordinate
  for (const auto l : lhsList) {
    if (lhsList.size() > 1) {
      gradients.row(l) -= (gradients.row(l).cwiseProduct(reactionCoordinate.row(l)).sum()) * reactionCoordinate.row(l);
    }
    else {
      gradients.row(l).setZero();
    }
    // replace gradient along reaction coordinate with fixed force only for lhs
    if (addForce && this->movableSide == "lhs") {
      double factor = (attractive) ? 1.0 : -1.0;
      gradients.row(l) += factor * totalForceNorm * reactionCoordinate.row(l);
    }
  }
  for (const auto r : rhsList) {
    if (rhsList.size() > 1) {
      gradients.row(r) -= (gradients.row(r).cwiseProduct(reactionCoordinate.row(r)).sum()) * reactionCoordinate.row(r);
    }
    else {
      gradients.row(r).setZero();
    }
    // replace gradient along reaction coordinate with fixed force only for rhs
    if (addForce && this->movableSide == "rhs") {
      double factor = (attractive) ? 1.0 : -1.0;
      gradients.row(r) += factor * totalForceNorm * reactionCoordinate.row(r);
    }
  }
  // replace gradient along reaction coordinate with fixed force for both sides
  if (addForce && this->movableSide == "both") {
    double factor = (attractive) ? 0.5 : -0.5;
    gradients += factor * totalForceNorm * reactionCoordinate;
  }
  // Apply Cartesian constraints
  if (!this->fixedAtoms.empty()) {
    for (const auto& a : fixedAtoms) {
      gradients.row(a).setZero();
    }
  }
}

Displacement NtOptimizer::centerToCenterVector(const PositionCollection& positions) const {
  Position lhsCenter = Position::Zero();
  Position rhsCenter = Position::Zero();
  for (const auto l : lhsList) {
    lhsCenter += positions.row(l);
  }
  lhsCenter.array() /= static_cast<double>(lhsList.size());
  for (const auto r : rhsList) {
    rhsCenter += positions.row(r);
  }
  rhsCenter.array() /= static_cast<double>(rhsList.size());
  Displacement c2c = rhsCenter - lhsCenter;
  return c2c;
}

bool NtOptimizer::convergedOptimization(const AtomCollection& atoms) const {
  const auto& positions = atoms.getPositions();
  double c2cNorm = this->centerToCenterVector(positions).norm();
  if (this->attractive) {
    if (c2cNorm < this->check.attractiveStop) {
      return true;
    }
    for (const auto& l : this->lhsList) {
      for (const auto& r : this->rhsList) {
        Eigen::Vector3d dir = positions.row(l) - positions.row(r);
        const double dist = dir.norm();
        const double r12cov =
            ElementInfo::covalentRadius(atoms.getElement(l)) + ElementInfo::covalentRadius(atoms.getElement(r));
        if (dist < this->check.attractiveStop * r12cov) {
          return true;
        }
      }
    }
    return false;
  }
  bool done = true;
  for (const auto& l : this->lhsList) {
    for (const auto& r : this->rhsList) {
      Eigen::Vector3d dir = positions.row(l) - positions.row(r);
      const double dist = dir.norm();
      const double r12cov =
          ElementInfo::covalentRadius(atoms.getElement(l)) + ElementInfo::covalentRadius(atoms.getElement(r));
      if (dist < this->check.repulsiveStop * r12cov) {
        done = false;
        break;
      }
    }
  }
  if (c2cNorm <= this->check.repulsiveStop) {
    done = false;
  }
  return done;
}

PositionCollection NtOptimizer::extractTsGuess() const {
  // Apply a Savitzky-Golay filter (width 5) and take the 1st derivative
  std::vector<double> filteredValues(_values);
  std::vector<double> filteredGradients(_values.size());
  for (int i = 0; i < this->filterPasses; i++) {
    // prepend and post pend data
    std::vector<double> tmp;
    tmp.reserve(filteredValues.size() + 4);
    tmp.push_back(filteredValues[0]);
    tmp.push_back(filteredValues[0]);
    tmp.insert(tmp.begin() + 2, filteredValues.begin(), filteredValues.end());
    tmp.push_back(filteredValues[filteredValues.size() - 1]);
    tmp.push_back(filteredValues[filteredValues.size() - 1]);
    for (unsigned int j = 2; j < filteredValues.size() + 2; j++) {
      const double a0 = (-3.0 * tmp[j - 2] + 12.0 * tmp[j - 1] + 17.0 * tmp[j] + 12.0 * tmp[j + 1] - 3.0 * tmp[j + 2]) / 35.0;
      const double a1 = (tmp[j - 2] - 8.0 * tmp[j - 1] + 8.0 * tmp[j + 1] - 1.0 * tmp[j + 2]) / 12.0;
      filteredValues[j - 2] = a0;
      filteredGradients[j - 2] = a1;
    }
  }

  // Find all maxima of energy curve as indicated by a plus to minus zero pass (read from the left)
  // Search from the back if the force is attractive
  std::vector<int> maximaList;
  if (this->attractive) {
    for (int i = static_cast<int>(filteredGradients.size()) - 2; i > 0; i--) {
      if (filteredGradients[i] >= 0.0 && filteredGradients[i + 1] < 0.0) {
        maximaList.push_back((fabs(filteredGradients[i]) < fabs(filteredGradients[i + 1])) ? i : i + 1);
      }
    }
  }
  else {
    const int N = static_cast<int>(_values.size()) - 1;
    for (int i = 0; i < N; i++) {
      if (filteredGradients[i + 1] <= 0.0 && filteredGradients[i] > 0.0) {
        maximaList.push_back((fabs(filteredGradients[i]) < fabs(filteredGradients[i + 1])) ? i : i + 1);
      }
    }
  }

  if (maximaList.empty()) {
    throw std::runtime_error("No transition state guess was found in Newton Trajectory scan.");
  }
  // Extract TS guess from point with highest energy
  double maxValue = -std::numeric_limits<double>::max();
  int maxIndex = -1;
  for (auto& i : maximaList) {
    if (_values[i] > maxValue) {
      maxIndex = i;
      maxValue = _values[i];
    }
  }

  return _trajectory[maxIndex];
}

void NtOptimizer::updateCoordinates(PositionCollection& coordinates, const AtomCollection& atoms,
                                    const GradientCollection& gradients) const {
  // transform coordinates and gradients if transformation was specified
  if (this->coordinateSystem == CoordinateSystem::Internal) {
    try {
      auto transformation = std::make_shared<InternalCoordinates>(atoms);
      auto internalCoordinates = transformation->coordinatesToInternal(coordinates);
      auto internalGradients = transformation->gradientsToInternal(gradients);
      internalCoordinates -= optimizer.factor * internalGradients;
      coordinates = transformation->coordinatesToCartesian(internalCoordinates);
    }
    catch (const InternalCoordinatesException& e) {
      // if true internals fail, fall back to Cartesians without rotation and translation
      auto transformation = std::make_shared<InternalCoordinates>(atoms, true);
      auto internalCoordinates = transformation->coordinatesToInternal(coordinates);
      auto internalGradients = transformation->gradientsToInternal(gradients);
      internalCoordinates -= optimizer.factor * internalGradients;
      coordinates = transformation->coordinatesToCartesian(internalCoordinates);
    }
  }
  else if (this->coordinateSystem == CoordinateSystem::CartesianWithoutRotTrans) {
    auto transformation = std::make_shared<InternalCoordinates>(atoms, true);
    auto internalCoordinates = transformation->coordinatesToInternal(coordinates);
    auto internalGradients = transformation->gradientsToInternal(gradients);
    internalCoordinates -= optimizer.factor * internalGradients;
    coordinates = transformation->coordinatesToCartesian(internalCoordinates);
  }
  else if (this->coordinateSystem == CoordinateSystem::Cartesian) {
    coordinates -= optimizer.factor * gradients;
  }
  else {
    throw std::runtime_error("Unknown coordinate system, please check your '" + std::string(ntCoordinateSystemKey) + "' input.");
  }
}

void NtOptimizer::addSettingsDescriptors(UniversalSettings::DescriptorCollection& /* collection */) const {
  // Not implemented yet
  throw std::runtime_error("You reached a function in the NTOptimizer that should not be called.");
}

void NtOptimizer::setSettings(const Settings& settings) {
  // optimizer.applySettings(settings);
  if (!settings.valid()) {
    settings.throwIncorrectSettings();
  }
  this->optimizer.factor = settings.getDouble(NtOptimizer::ntSdFactorKey);
  this->check.maxIter = settings.getInt(NtOptimizer::ntMaxIterKey);
  this->check.repulsiveStop = settings.getDouble(NtOptimizer::ntRepulsiveStopKey);
  this->check.attractiveStop = settings.getDouble(NtOptimizer::ntAttractiveStopKey);
  this->rhsList = settings.getIntList(NtOptimizer::ntRHSListKey);
  this->lhsList = settings.getIntList(NtOptimizer::ntLHSListKey);
  this->attractive = settings.getBool(NtOptimizer::ntAttractiveKey);
  this->totalForceNorm = settings.getDouble(NtOptimizer::ntTotalForceNormKey);
  this->coordinateSystem =
      CoordinateSystemInterpreter::getCoordinateSystemFromString(settings.getString(NtOptimizer::ntCoordinateSystemKey));
  this->useMicroCycles = settings.getBool(NtOptimizer::ntUseMicroCycles);
  this->fixedNumberOfMicroCycles = settings.getBool(NtOptimizer::ntFixedNumberOfMicroCycles);
  this->numberOfMicroCycles = settings.getInt(NtOptimizer::ntNumberOfMicroCycles);
  this->filterPasses = settings.getInt(NtOptimizer::ntFilterPasses);
  this->fixedAtoms = settings.getIntList(NtOptimizer::ntFixedAtomsKey);
  this->movableSide = settings.getString(NtOptimizer::ntMovableSide);

  // Check whether constraints and coordinate transformations are both switched on:
  if (!this->fixedAtoms.empty() && this->coordinateSystem != CoordinateSystem::Cartesian) {
    throw std::logic_error("Cartesian constraints cannot be set when using coordinate transformations! Set "
                           "'nt_coordinate_system' to 'cartesian'.");
  }
}

void NtOptimizer::applySettings(const Settings& settings) {
  this->setSettings(settings);
}

Settings NtOptimizer::getSettings() const {
  return NtOptimizerSettings(*this);
}

} // namespace Utils
} // namespace Scine
