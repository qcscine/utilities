/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MolecularTrajectory.h"

namespace Scine {
namespace Utils {

MolecularTrajectory::MolecularTrajectory(const ElementTypeCollection& elements) {
  elements_ = elements;
}

MolecularTrajectory::MolecularTrajectory(double minimumRmsdForAddition)
  : MolecularTrajectory(elements_, minimumRmsdForAddition) {
}

MolecularTrajectory::MolecularTrajectory(const ElementTypeCollection& elements, double minimumRmsdForAddition) {
  elements_ = elements;
  minMeanSquareDeviation_ = minimumRmsdForAddition * minimumRmsdForAddition;
  respectMinRmsd_ = minimumRmsdForAddition > 1e-12;
}

void MolecularTrajectory::setElementType(int index, ElementType e) {
  elements_[index] = e;
}

void MolecularTrajectory::setElementTypes(const ElementTypeCollection& elements) {
  assert(resettingElementTypeCollectionIsAllowed(elements));
  elements_ = elements;
}

const ElementTypeCollection& MolecularTrajectory::getElementTypes() const {
  return elements_;
}

void MolecularTrajectory::setEnergies(const EnergyContainer& energies) {
  if (static_cast<int>(energies.size()) != size()) {
    throw std::runtime_error(
        "The number of energies does not match the number of structures in this molecular trajectory.");
  }
  energies_ = energies;
}

MolecularTrajectory::EnergyContainer MolecularTrajectory::getEnergies() const {
  return energies_;
}

std::vector<PeriodicBoundaries> MolecularTrajectory::getPbcs() const {
  std::vector<PeriodicBoundaries> pbcs;
  std::transform(pbcs_.begin(), pbcs_.end(), std::back_inserter(pbcs),
                 [&](Eigen::Matrix3d m) -> PeriodicBoundaries { return PeriodicBoundaries(m); });
  return pbcs;
}

void MolecularTrajectory::setPbcs(const PbcContainer& pbcs) {
  if (static_cast<int>(pbcs.size()) != size()) {
    throw std::runtime_error(
        "The number of pbcs does not match the number of structures in this molecular trajectory.");
  }
  pbcs_ = pbcs;
}

void MolecularTrajectory::setPbcs(const std::vector<PeriodicBoundaries>& pbcs) {
  PbcContainer newPbcs;
  std::transform(pbcs.begin(), pbcs.end(), std::back_inserter(newPbcs),
                 [&](PeriodicBoundaries pbc) -> Eigen::Matrix3d { return pbc.getCellMatrix(); });
  setPbcs(newPbcs);
}

void MolecularTrajectory::clearEnergies() {
  energies_.clear();
}

void MolecularTrajectory::clearPbcs() {
  pbcs_.clear();
}

void MolecularTrajectory::clear() {
  structureVector_.clear();
  energies_.clear();
  pbcs_.clear();
}

void MolecularTrajectory::resize(int n) {
  structureVector_.resize(static_cast<Container::size_type>(n));
  energies_.resize(static_cast<Container::size_type>(n));
  pbcs_.resize(static_cast<Container::size_type>(n));
}

void MolecularTrajectory::push_back(PositionCollection p) {
  assert(additionOfPositionCollectionIsAllowed(p));
  if (!energies_.empty()) {
    throw std::runtime_error("Energy container is not empty. Clear energies first, before a structure only is added.");
  }
  if (!pbcs_.empty()) {
    throw std::runtime_error("Pbc container is not empty. Clear pbcs first, before a structure only is added.");
  }
  if (additionOfPositionCollectionIsAllowedBasedOnRmsd(p)) {
    structureVector_.push_back(std::move(p));
  }
}

void MolecularTrajectory::push_back(PositionCollection p, double e) {
  assert(additionOfPositionCollectionIsAllowed(p));
  if (static_cast<int>(energies_.size()) != size()) {
    throw std::runtime_error(
        "The number of energies does not match the number of structures in this molecular trajectory.");
  }
  if (additionOfPositionCollectionIsAllowedBasedOnRmsd(p)) {
    structureVector_.push_back(std::move(p));
    energies_.push_back(e);
  }
}

void MolecularTrajectory::push_back(PositionCollection p, PeriodicBoundaries pbc) {
  assert(additionOfPositionCollectionIsAllowed(p));
  if (static_cast<int>(pbcs_.size()) != size()) {
    throw std::runtime_error(
        "The number of pbcs does not match the number of structures in this molecular trajectory.");
  }
  if (additionOfPositionCollectionIsAllowedBasedOnRmsd(p)) {
    structureVector_.push_back(std::move(p));
    pbcs_.push_back(pbc.getCellMatrix());
  }
}

void MolecularTrajectory::push_back(PositionCollection p, double e, PeriodicBoundaries pbc) {
  assert(additionOfPositionCollectionIsAllowed(p));
  auto s = size();
  if (static_cast<int>(energies_.size()) != s) {
    throw std::runtime_error(
        "The number of energies does not match the number of structures in this molecular trajectory.");
  }
  if (static_cast<int>(pbcs_.size()) != s) {
    throw std::runtime_error(
        "The number of pbcs does not match the number of structures in this molecular trajectory.");
  }
  if (additionOfPositionCollectionIsAllowedBasedOnRmsd(p)) {
    structureVector_.push_back(std::move(p));
    energies_.push_back(e);
    pbcs_.push_back(pbc.getCellMatrix());
  }
}

bool MolecularTrajectory::empty() const {
  return structureVector_.empty();
}

int MolecularTrajectory::size() const {
  return static_cast<int>(structureVector_.size());
}

int MolecularTrajectory::molecularSize() const {
  return elements_.size();
}

MolecularTrajectory::iterator MolecularTrajectory::begin() {
  return structureVector_.begin();
}

MolecularTrajectory::const_iterator MolecularTrajectory::begin() const {
  return structureVector_.begin();
}

MolecularTrajectory::iterator MolecularTrajectory::erase(iterator position) {
  return structureVector_.erase(position);
}

MolecularTrajectory::iterator MolecularTrajectory::end() {
  return structureVector_.end();
}

MolecularTrajectory::const_iterator MolecularTrajectory::end() const {
  return structureVector_.end();
}

MolecularTrajectory::reference MolecularTrajectory::operator[](int i) {
  return structureVector_[i];
}

MolecularTrajectory::const_reference MolecularTrajectory::operator[](int i) const {
  return structureVector_[i];
}

MolecularTrajectory::reference MolecularTrajectory::at(int i) {
  return structureVector_.at(static_cast<Container::size_type>(i));
}

MolecularTrajectory::const_reference MolecularTrajectory::at(int i) const {
  return structureVector_.at(static_cast<Container::size_type>(i));
}

MolecularTrajectory::reference MolecularTrajectory::front() {
  return structureVector_.front();
}

MolecularTrajectory::const_reference MolecularTrajectory::front() const {
  return structureVector_.front();
}

MolecularTrajectory::reference MolecularTrajectory::back() {
  return structureVector_.back();
}

MolecularTrajectory::const_reference MolecularTrajectory::back() const {
  return structureVector_.back();
}

const MolecularTrajectory& MolecularTrajectory::operator*=(double f) {
  for (auto& s : structureVector_) {
    s *= f;
  }
  for (auto& pbc : pbcs_) {
    pbc *= f;
  }
  return *this;
}

const MolecularTrajectory& MolecularTrajectory::operator/=(double f) {
  for (auto& s : structureVector_) {
    s /= f;
  }
  for (auto& pbc : pbcs_) {
    pbc /= f;
  }
  return *this;
}

MolecularTrajectory MolecularTrajectory::operator*(double f) const {
  MolecularTrajectory t = *this;
  t *= f;
  return t;
}

MolecularTrajectory MolecularTrajectory::operator/(double f) const {
  MolecularTrajectory t = *this;
  t /= f;
  return t;
}

bool MolecularTrajectory::resettingElementTypeCollectionIsAllowed(const ElementTypeCollection& ec) const {
  bool hasSameSizeAsPreviousElementTypeCollection = static_cast<int>(ec.size()) == molecularSize();
  bool noPositionsArePresent = empty();
  bool hasSameSizeAsPresentPositions = empty() || (structureVector_.front().rows() == static_cast<Eigen::Index>(ec.size()));
  return hasSameSizeAsPreviousElementTypeCollection || noPositionsArePresent || hasSameSizeAsPresentPositions;
}

bool MolecularTrajectory::additionOfPositionCollectionIsAllowed(const PositionCollection& p) const {
  bool validWithRespectToOtherPositionCollections = empty() || p.size() == structureVector_.front().size();
  bool validWithRespectToElementTypeCollection = p.rows() == molecularSize();

  bool validGivenThatElementCollectionIsAlreadySet = validWithRespectToElementTypeCollection;
  bool validGivenThatElementCollectionIsNotSet = (molecularSize() == 0 && validWithRespectToOtherPositionCollections);

  return validGivenThatElementCollectionIsAlreadySet || validGivenThatElementCollectionIsNotSet;
}

bool MolecularTrajectory::additionOfPositionCollectionIsAllowedBasedOnRmsd(const PositionCollection& p) const {
  if (!respectMinRmsd_ || structureVector_.empty()) {
    return true;
  }
  PositionCollection lastPos = structureVector_.back();
  double meanSquareDeviation = (lastPos - p).rowwise().squaredNorm().sum() / lastPos.rows();
  return meanSquareDeviation > minMeanSquareDeviation_;
}

void MolecularTrajectory::setResidues(const ResidueCollection& residues) {
  if (residues.size() != unsigned(this->molecularSize())) {
    throw std::runtime_error(
        "The number of residue information must be identical to the number of atoms in the collection");
  }
  residues_ = residues;
}

const ResidueCollection& MolecularTrajectory::getResidues() const {
  return residues_;
}

} /* namespace Utils */
} /* namespace Scine */
