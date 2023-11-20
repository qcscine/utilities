/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/GeometryOptimization/NtOptimizer2.h"
#include "Utils/CalculatorBasics/CalculationRoutines.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/GeometryOptimization/NtOptimizer2Settings.h"
#include "Utils/Optimizer/GradientBased/Bfgs.h"
#include <array>
#include <boost/exception/diagnostic_information.hpp>
#include <cmath>
#include <valarray>

namespace Scine {
namespace Utils {

int NtOptimizer2::optimize(AtomCollection& atoms, Core::Log& log) {
  this->sanityCheck(atoms);
  this->setReactiveAtomsList();
  this->setConstraintsMap(atoms);
  // Configure Calculator
  _calculator.setStructure(atoms);
  _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::BondOrderMatrix);
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
        _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::BondOrderMatrix);
        atoms.setPositions(coordinates);
        Results results = CalculationRoutines::calculateWithCatch(_calculator, log, "Calculation in NT optimization failed.");
        value = results.get<Property::Energy>();
        // Apply Cartesian constraints
        auto bos = results.get<Property::BondOrderMatrix>();
        auto gradientMatrix = results.get<Property::Gradients>();
        this->updateGradients(atoms, value, gradientMatrix, bos, cycle);
        gradients = Eigen::Map<const Eigen::VectorXd>(gradientMatrix.data(), nAtoms * 3);
      };
      // Run micro iterations
      Eigen::VectorXd positions = Eigen::Map<const Eigen::VectorXd>(atoms.getPositions().data(), atoms.size() * 3);
      try {
        bfgs.optimize(positions, microIterUpdate, microIterCheck, log);
      }
      catch (...) {
        if (cycle > 5) {
          coordinates = this->extractTsGuess();
          atoms.setPositions(coordinates);
          _calculator.modifyPositions(coordinates);
          return cycle;
        }
        throw std::runtime_error("NT micro iterations failed:\n" + boost::current_exception_diagnostic_information());
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
    auto bos = results.get<Property::BondOrderMatrix>();
    this->triggerObservers(cycle, value, Eigen::Map<const Eigen::VectorXd>(coordinates.data(), nAtoms * 3));
    // Evaluate additional force
    this->updateGradients(atoms, value, gradients, bos, cycle, true);
    // Add energy to energy vector
    _values.push_back(value);
    // Add coordinates to trajectory
    _trajectory.push_back(coordinates);
    // Check convergence
    if (this->convergedOptimization(atoms, bos)) {
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
  // If the maximum number of iterations is reached try to extract a guess
  coordinates = this->extractTsGuess();
  atoms.setPositions(coordinates);
  _calculator.modifyPositions(coordinates);
  return cycle;
}

void NtOptimizer2::sanityCheck(const AtomCollection& atoms) const {
  if (associationList.empty() && dissociationList.empty()) {
    throw std::logic_error("Called NT optimization without atoms to apply a force to.");
  }
  // check if both list have even numbers of atom indices
  if (associationList.size() % 2 != 0) {
    throw std::logic_error("Called NT optimization with uneven number of atom indices for bonds.");
  }
  if (dissociationList.size() % 2 != 0) {
    throw std::logic_error("Called NT optimization with uneven number of atom indices for bonds.");
  }
  // check for negative indices
  if (std::any_of(associationList.begin(), associationList.end(), [&](int i) { return i < 0; })) {
    throw std::logic_error(
        "At least one of the given indices in the 'nt_associations' is negative, which is not possible.");
  }
  if (std::any_of(dissociationList.begin(), dissociationList.end(), [&](int i) { return i < 0; })) {
    throw std::logic_error(
        "At least one of the given indices in the 'nt_dissociations' is negative, which is not possible.");
  }
  // check for too large indices
  const int nAtoms = atoms.size();
  for (const auto i : associationList) {
    if (i >= nAtoms) {
      throw std::logic_error("Index '" + std::to_string(i) + "' in 'nt_associations' too large for " +
                             std::to_string(nAtoms) + " atoms.");
    }
  }
  for (const auto i : dissociationList) {
    if (i >= nAtoms) {
      throw std::logic_error("Index '" + std::to_string(i) + "' in 'nt_dissociations' too large for " +
                             std::to_string(nAtoms) + " atoms.");
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
      movableReactionIndicesList = associationList;
      movableReactionIndicesList.insert(movableReactionIndicesList.end(), dissociationList.begin(), dissociationList.end());
      if (std::find(movableReactionIndicesList.begin(), movableReactionIndicesList.end(), a) !=
          movableReactionIndicesList.end()) {
        throw std::logic_error("The atom index " + std::to_string(a) +
                               " was specified to be constrained although "
                               "it was declared as part of the reaction coordinate");
      }
    }
  }
  if (std::find(possibleExtractionOptions.begin(), possibleExtractionOptions.end(), extractionCriterion) ==
      possibleExtractionOptions.end()) {
    std::string options = "\nValid options are:\n";
    for (const auto& criterion : possibleExtractionOptions) {
      options += "'" + criterion + "'\n";
    }
    throw std::logic_error("Extract criterion " + extractionCriterion + " is not a valid option." + options);
  }
}

void NtOptimizer2::updateGradients(const AtomCollection& atoms, const double& /* energy */, GradientCollection& gradients,
                                   const BondOrderCollection& bondOrders, int cycle, bool addForce) {
  const auto& positions = atoms.getPositions();
  eliminateReactiveAtomsGradients(positions, gradients);
  auto reactions = inferReactions(bondOrders);
  ReactionMapping associations = reactions.first;
  ReactionMapping dissociations = reactions.second;
  // define reaction coordinate
  GradientCollection reactionCoordinate(gradients);
  reactionCoordinate.setZero();
  double maxScale = 0.0;
  // association
  for (const auto& [left, right] : associations) {
    double r12cov = NtUtils::smallestCovalentRadius(atoms, left) + NtUtils::smallestCovalentRadius(atoms, right);
    Displacement c2c = NtUtils::centerToCenterVector(positions, left, right);
    double dist = c2c.norm();
    double bo = 0.0;
    for (const auto& l : left) {
      for (const auto& r : right) {
        bo += bondOrders.getOrder(l, r);
      }
    }
    /* Check if one should keep pushing along the reaction coordinate */
    /* Keep pushing until 10% above attractive bond order stop OR distance 10% below distance stop */
    if (bo > 1.1 * check.attractiveBondOrderStop || dist < 0.9 * this->check.attractiveDistanceStop * r12cov) {
      if (_firstCoordinateReachedIndex == -1) {
        _firstCoordinateReachedIndex = cycle;
      }
      continue;
    }
    double scale = dist - r12cov;
    if (scale < 0) {
      scale = 0.1;
    }
    maxScale = std::max(maxScale, scale);
    c2c *= scale / dist;
    for (const auto& l : left) {
      for (const auto& r : right) {
        reactionCoordinate.row(l) += c2c;
        reactionCoordinate.row(r) -= c2c;
      }
    }
  }
  // dissociation
  for (const auto& [left, right] : dissociations) {
    double bo = 0.0;
    Displacement c2c = NtUtils::centerToCenterVector(positions, left, right);
    double dist = c2c.norm();
    for (const auto& l : left) {
      for (const auto& r : right) {
        bo += bondOrders.getOrder(l, r);
      }
    }
    /* Check if one should keep pulling along the reaction coordinate */
    /* Keep pulling until 30% bellow repulsive bond order stop */
    if (bo < 0.7 * this->check.repulsiveBondOrderStop) {
      if (_firstCoordinateReachedIndex == -1) {
        _firstCoordinateReachedIndex = cycle;
      }
      continue;
    }
    const double scale = 0.8;
    maxScale = std::max(maxScale, scale);
    c2c *= scale / dist;
    for (const auto& l : left) {
      for (const auto& r : right) {
        reactionCoordinate.row(l) -= c2c;
        reactionCoordinate.row(r) += c2c;
      }
    }
  }
  if (maxScale != 0.0) {
    reactionCoordinate *= 0.5 * (totalForceNorm / maxScale);
  }

  /* Update gradients */
  // replace gradient along reaction coordinate with fixed force for both sides
  if (addForce) {
    gradients += reactionCoordinate;
  }
  // Apply Cartesian constraints
  if (!this->fixedAtoms.empty()) {
    for (const auto& a : fixedAtoms) {
      gradients.row(a).setZero();
    }
  }
}

void NtOptimizer2::eliminateReactiveAtomsGradients(const PositionCollection& positions, GradientCollection& gradients) const {
  // remove gradient along reaction coordinates
  /* Loop over all reactive atoms */
  for (auto const& reactiveAtomIndex : _reactiveAtomsList) {
    /* Build a normalized matrix of the constraints */
    PositionCollection constraintsMatrix;
    for (unsigned int j = 0; j < _constraintsMap[reactiveAtomIndex].size() / 2; j++) {
      const auto& l = _constraintsMap[reactiveAtomIndex][2 * j];
      const auto& r = _constraintsMap[reactiveAtomIndex][2 * j + 1];
      constraintsMatrix.conservativeResize(constraintsMatrix.rows() + 1, 3);
      constraintsMatrix.row(constraintsMatrix.rows() - 1) =
          (positions.row(l) - positions.row(r)) / (positions.row(l) - positions.row(r)).norm();
    }
    /* Solve square axes matrix with SelfAdjointEigenSolver */
    Eigen::MatrixXd constraintsMatrixSquare = constraintsMatrix * constraintsMatrix.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes;
    saes.compute(constraintsMatrixSquare);
    Eigen::VectorXd saesEigenVal = saes.eigenvalues();
    Eigen::MatrixXd saesEigenVec = saes.eigenvectors();
    PositionCollection constraintsMatrixOrtho = (constraintsMatrix.transpose() * saesEigenVec).transpose();
    /* Determine rank of matrix (eigenvalues > 1e-6) */
    unsigned int rank = 0;
    long lastRowIndex = constraintsMatrixOrtho.rows() - 1;
    /* Determine rank of matrix */
    for (long i = constraintsMatrixOrtho.rows() - 1; i > -1; i--) {
      /* Break loop over eigenvalues if it is < 1e-6 */
      if (saesEigenVal.col(0)[i] < 1e-6) {
        break;
      }
      rank += 1;
    }
    /* Manipulate gradients of reactive atom */
    /* Can not be orthogonal to more than 2 axis, hence set the gradient to 0 */
    if (rank > 2) {
      gradients.row(reactiveAtomIndex).setZero();
    }
    /* Project gradient on the cross product of 2 linear independent axes */
    else if (rank == 2) {
      Eigen::Vector3d a2a0 = constraintsMatrixOrtho.row(lastRowIndex) / (constraintsMatrix.row(lastRowIndex).norm());
      Eigen::Vector3d a2a1 = constraintsMatrixOrtho.row(lastRowIndex - 1) / (constraintsMatrix.row(lastRowIndex - 1).norm());
      Eigen::Vector3d proj = a2a0.cross(a2a1);
      proj /= proj.norm();
      gradients.row(reactiveAtomIndex) = (gradients.row(reactiveAtomIndex).dot(proj)) * proj;
    }
    /* Subtract projection along the axis from gradient */
    else if (rank == 1) {
      Eigen::Vector3d a2a = constraintsMatrixOrtho.row(lastRowIndex) / (constraintsMatrixOrtho.row(lastRowIndex).norm());
      gradients.row(reactiveAtomIndex) -= (gradients.row(reactiveAtomIndex).dot(a2a)) * a2a;
    }
  }
}

std::pair<NtOptimizer2::ReactionMapping, NtOptimizer2::ReactionMapping>
NtOptimizer2::inferReactions(const BondOrderCollection& bondOrders) const {
  std::vector<ReactionMapping> maps;
  // we do the exact same thing for associations and dissociations
  std::vector<std::vector<int>> reactionLists = {associationList, dissociationList};
  for (const auto& reactionList : reactionLists) {
    // we restructure the single list with included logic into list of pairs
    std::vector<std::pair<int, int>> reactions;
    for (unsigned int i = 0; i < reactionList.size() / 2; i++) {
      const auto& left = reactionList[2 * i];
      const auto& right = reactionList[2 * i + 1];
      reactions.emplace_back(left, right);
    }
    // now we need to combine duplicates
    // we pretty much build now a map that allows duplicate keys
    // reason is, if we have a duplicate index X on the left side, we want to distinguish if X takes part in genuinely
    // different reaction coordinates, e.g. X - Y and X - Z, or if we have an eta bond reaction coordinate X - Y-Z
    // we distinguish this in our data structure with:
    // <X, {Y}>
    // <X, {Z}>
    // for different reaction coordinates and
    // <X, {Y, Z}>
    // for an eta bond reaction coordinate. The value can have any length and is not limited to two
    // we determine the case based on the fact if Y and Z are part of the same molecule (second case) or not (first case)
    std::vector<std::pair<int, std::vector<int>>> rightCorrectedReactions;
    std::vector<int> sanitizedLeftIndices;
    for (unsigned long i = 0; i < reactions.size(); ++i) {
      if (std::find(sanitizedLeftIndices.begin(), sanitizedLeftIndices.end(), reactions[i].first) !=
          sanitizedLeftIndices.end()) {
        continue; // no need to check same index again
      }
      sanitizedLeftIndices.push_back(reactions[i].first);
      std::vector<int> values = {reactions[i].second};
      for (unsigned long j = i + 1; j < reactions.size(); ++j) {
        if (reactions[i].first == reactions[j].first) {
          values.push_back(reactions[j].second);
        }
      }
      if (values.size() > 1) {
        // we had a duplicate, let's check if the shared values are part of the same molecule
        std::vector<std::vector<int>> molecules = NtUtils::connectedNuclei(values, bondOrders);
        // add entry for each molecule, since each represents a genuine reaction coordinate
        for (const auto& molecule : molecules) {
          rightCorrectedReactions.emplace_back(reactions[i].first, molecule);
        }
      }
      else {
        // no duplicate simply copy into new data structure
        rightCorrectedReactions.emplace_back(reactions[i].first, values);
      }
    }
    // now we need to do the exact same thing for the right side and check if there are indices on the right side
    // or lists of indices that have different reaction partners and check if they are again different reaction
    // coordinates or eta bond formations
    ReactionMapping allCorrectedReactions;
    std::vector<std::vector<int>> sanitizedRightIndices;
    for (unsigned long i = 0; i < rightCorrectedReactions.size(); ++i) {
      if (std::find(sanitizedRightIndices.begin(), sanitizedRightIndices.end(), rightCorrectedReactions[i].second) !=
          sanitizedRightIndices.end()) {
        continue;
      }
      sanitizedRightIndices.push_back(rightCorrectedReactions[i].second);
      std::vector<int> values = {rightCorrectedReactions[i].first};
      for (unsigned long j = i + 1; j < rightCorrectedReactions.size(); ++j) {
        if (rightCorrectedReactions[i].second == rightCorrectedReactions[j].second) {
          values.push_back(rightCorrectedReactions[j].first);
        }
      }
      if (values.size() > 1) {
        std::vector<std::vector<int>> molecules = NtUtils::connectedNuclei(values, bondOrders);
        for (const auto& molecule : molecules) {
          allCorrectedReactions.emplace_back(molecule, rightCorrectedReactions[i].second);
        }
      }
      else {
        allCorrectedReactions.emplace_back(values, rightCorrectedReactions[i].second);
      }
    }
    maps.push_back(allCorrectedReactions);
  }
  return std::make_pair(maps[0], maps[1]);
}

bool NtOptimizer2::convergedOptimization(const AtomCollection& atoms, const BondOrderCollection& bondOrders) const {
  auto reactions = inferReactions(bondOrders);
  ReactionMapping associations = reactions.first;
  ReactionMapping dissociations = reactions.second;
  const auto& positions = atoms.getPositions();
  // associations
  for (const auto& [left, right] : associations) {
    double r12cov = NtUtils::smallestCovalentRadius(atoms, left) + NtUtils::smallestCovalentRadius(atoms, right);
    Displacement c2c = NtUtils::centerToCenterVector(positions, left, right);
    double dist = c2c.norm();
    double bo = 0.0;
    for (const auto& l : left) {
      for (const auto& r : right) {
        bo += bondOrders.getOrder(l, r);
      }
    }
    if (bo < this->check.attractiveBondOrderStop && dist > this->check.attractiveDistanceStop * r12cov) {
      return false;
    }
  }
  // dissociations
  for (const auto& [left, right] : dissociations) {
    double bo = 0.0;
    for (const auto& l : left) {
      for (const auto& r : right) {
        bo += bondOrders.getOrder(l, r);
      }
    }
    if (bo > this->check.repulsiveBondOrderStop) {
      return false;
    }
  }
  return true;
}

PositionCollection NtOptimizer2::extractTsGuess() const {
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
  for (int i = static_cast<int>(filteredGradients.size()) - 2; i > 0; i--) {
    if (filteredGradients[i] >= 0.0 && filteredGradients[i + 1] < 0.0) {
      maximaList.push_back((fabs(filteredGradients[i]) < fabs(filteredGradients[i + 1])) ? i : i + 1);
    }
  }

  if (maximaList.empty()) {
    throw std::runtime_error("No transition state guess was found in Newton Trajectory scan.");
  }
  // Extract TS guess according to given criterion
  if (extractionCriterion == SettingsNames::Optimizations::Nt2::extractFirst) {
    return _trajectory[maximaList.back()]; // back because maximalist starts from back
  }
  // get the highest value if wished or coordinate was never stopped
  if (extractionCriterion == SettingsNames::Optimizations::Nt2::extractHighest || _firstCoordinateReachedIndex == -1) {
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
  int firstMaximumAfterCoordinateReached = -1;
  // maximaList lists maxima from back
  for (auto& maximum : maximaList) {
    if (maximum < _firstCoordinateReachedIndex) {
      return _trajectory[maximum];
    }
    firstMaximumAfterCoordinateReached = maximum;
  }
  assert(firstMaximumAfterCoordinateReached > -1);
  // got no maximum before coordinate was ended, so we pick the first after the coordinate was reached
  return _trajectory[firstMaximumAfterCoordinateReached];
}

void NtOptimizer2::updateCoordinates(PositionCollection& coordinates, const AtomCollection& atoms,
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
    throw std::runtime_error("Unknown coordinate system, please check your '" +
                             std::string(SettingsNames::Optimizations::Nt2::coordinateSystem) + "' input.");
  }
}

void NtOptimizer2::addSettingsDescriptors(UniversalSettings::DescriptorCollection& /* collection */) const {
  // Not implemented yet
  throw std::runtime_error("You reached a function in the NTOptimizer that should not be called.");
}

void NtOptimizer2::setSettings(const Settings& settings) {
  if (!settings.valid()) {
    settings.throwIncorrectSettings();
  }
  this->optimizer.factor = settings.getDouble(SettingsNames::Optimizations::Nt2::sdFactor);
  this->check.maxIter = settings.getInt(SettingsNames::Optimizations::Nt2::maxIter);
  this->check.attractiveDistanceStop = settings.getDouble(SettingsNames::Optimizations::Nt2::attractiveStop);
  this->associationList = settings.getIntList(SettingsNames::Optimizations::Nt2::assList);
  this->dissociationList = settings.getIntList(SettingsNames::Optimizations::Nt2::dissList);
  this->totalForceNorm = settings.getDouble(SettingsNames::Optimizations::Nt2::totalForceNorm);
  this->coordinateSystem = CoordinateSystemInterpreter::getCoordinateSystemFromString(
      settings.getString(SettingsNames::Optimizations::Nt2::coordinateSystem));
  this->useMicroCycles = settings.getBool(SettingsNames::Optimizations::Nt2::useMicroCycles);
  this->fixedNumberOfMicroCycles = settings.getBool(SettingsNames::Optimizations::Nt2::fixedNumberOfMicroCycles);
  this->numberOfMicroCycles = settings.getInt(SettingsNames::Optimizations::Nt2::numberOfMicroCycles);
  this->filterPasses = settings.getInt(SettingsNames::Optimizations::Nt2::filterPasses);
  this->fixedAtoms = settings.getIntList(SettingsNames::Optimizations::Nt2::fixedAtoms);
  this->extractionCriterion = settings.getString(SettingsNames::Optimizations::Nt2::extractionCriterion);

  // Check whether constraints and coordinate transformations are both switched on:
  if (!this->fixedAtoms.empty() && this->coordinateSystem != CoordinateSystem::Cartesian) {
    throw std::logic_error("Cartesian constraints cannot be set when using coordinate transformations! Set "
                           "'nt_coordinate_system' to 'cartesian'.");
  }
}

void NtOptimizer2::applySettings(const Settings& settings) {
  this->setSettings(settings);
}

void NtOptimizer2::setReactiveAtomsList() {
  this->_reactiveAtomsList.clear();
  std::vector<int> l1 = this->associationList;
  std::vector<int> l2 = this->dissociationList;
  std::sort(l1.begin(), l1.end());
  std::sort(l2.begin(), l2.end());
  std::vector<int> r(l1.size() + l2.size());
  std::merge(l1.begin(), l1.end(), l2.begin(), l2.end(), std::back_inserter(this->_reactiveAtomsList));
  auto last = std::unique(this->_reactiveAtomsList.begin(), this->_reactiveAtomsList.end());
  this->_reactiveAtomsList.erase(last, this->_reactiveAtomsList.end());
}

void NtOptimizer2::setConstraintsMap(const AtomCollection& atoms) {
  this->_constraintsMap.resize(atoms.size());
  for (auto const& i : _reactiveAtomsList) {
    std::vector<int> matches = {};
    for (unsigned int j = 0; j < associationList.size(); j++) {
      if (associationList[j] == i) {
        if (j % 2 == 0) {
          matches.push_back(associationList[j]);
          matches.push_back(associationList[j + 1]);
        }
        else {
          matches.push_back(associationList[j - 1]);
          matches.push_back(associationList[j]);
        }
      }
    }
    for (unsigned int j = 0; j < dissociationList.size(); j++) {
      if (dissociationList[j] == i) {
        if (j % 2 == 0) {
          matches.push_back(dissociationList[j]);
          matches.push_back(dissociationList[j + 1]);
        }
        else {
          matches.push_back(dissociationList[j - 1]);
          matches.push_back(dissociationList[j]);
        }
      }
    }
    this->_constraintsMap[i] = matches;
  }
}

Settings NtOptimizer2::getSettings() const {
  return NtOptimizer2Settings(*this);
}

const std::vector<std::vector<int>>& NtOptimizer2::getConstraintsMap() {
  return _constraintsMap;
}

const std::vector<int>& NtOptimizer2::getReactiveAtomsList() {
  return _reactiveAtomsList;
}

namespace NtUtils {
std::vector<std::vector<int>> connectedNuclei(std::vector<int> indices, const BondOrderCollection& bondOrders) {
  std::sort(indices.begin(), indices.end()); // makes avoiding double counting easier
  std::vector<std::vector<int>> result;
  for (const auto& i : indices) {
    // check if i is already in one of the entries
    if (std::any_of(result.begin(), result.end(), [&](std::vector<int> bondedIndices) {
          return std::find(bondedIndices.begin(), bondedIndices.end(), i) != bondedIndices.end();
        })) {
      continue;
    }
    if (i == indices.back()) {
      // we reached last entry and this has never been included in the others, so it must be its own thing
      result.push_back({i});
      continue;
    }
    std::vector<int> connected;
    std::vector<int> nucleiToVisit = {i};
    while (!nucleiToVisit.empty()) {
      // take last value in list to crawl
      int current = nucleiToVisit.back();
      nucleiToVisit.pop_back();
      auto bondPartners = bondOrders.getBondPartners(current, 0.5);
      for (const auto& partner : bondPartners) {
        // add newly found indices to total list and list to further extend search
        if (std::find(connected.begin(), connected.end(), partner) == connected.end()) {
          connected.push_back(partner);
          nucleiToVisit.push_back(partner);
        }
      }
    }
    // now check for all others if they are within the list of bond partners i.e. same molecule
    std::vector<int> connectedIndices = {i};
    for (const auto& j : indices) {
      if (j > i && std::find(connected.begin(), connected.end(), j) != connected.end()) {
        connectedIndices.push_back(j);
      }
    }
    std::sort(connectedIndices.begin(), connectedIndices.end());
    result.push_back(connectedIndices);
  }
  return result;
}

Displacement centerToCenterVector(const PositionCollection& positions, const std::vector<int>& lhsList,
                                  const std::vector<int>& rhsList) {
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
  Displacement c2c = lhsCenter - rhsCenter;
  return c2c;
}

double smallestCovalentRadius(const AtomCollection& atoms, const std::vector<int>& indices) {
  if (indices.empty()) {
    throw std::runtime_error("Missing reactive atom indices");
  }
  double result = std::numeric_limits<double>::max();
  for (const auto& index : indices) {
    double rad = ElementInfo::covalentRadius(atoms.getElement(index));
    if (rad < result) {
      result = rad;
    }
  }
  return result;
}

} // namespace NtUtils
} // namespace Utils
} // namespace Scine
