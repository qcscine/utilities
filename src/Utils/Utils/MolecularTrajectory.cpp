/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MolecularTrajectory.h"

namespace Scine {
namespace Utils {

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
  if (energies.size() != size())
    throw std::runtime_error(
        "The number of energies does not match the number of structures in this molecular trajectory.");
  energies_ = energies;
}

MolecularTrajectory::EnergyContainer MolecularTrajectory::getEnergies() {
  return energies_;
}

void MolecularTrajectory::clearEnergies() {
  energies_.clear();
}

void MolecularTrajectory::clear() {
  structureVector_.clear();
  energies_.clear();
};

void MolecularTrajectory::resize(int n) {
  structureVector_.resize(static_cast<Container::size_type>(n));
  energies_.resize(static_cast<Container::size_type>(n));
}

void MolecularTrajectory::push_back(PositionCollection p) {
  assert(additionOfPositionCollectionIsAllowed(p));
  if (!energies_.empty())
    throw std::runtime_error("Energy container is not empty. Clear energies first, before a structure only is added.");
  structureVector_.push_back(std::move(p));
}

void MolecularTrajectory::push_back(PositionCollection p, double e) {
  assert(additionOfPositionCollectionIsAllowed(p));
  if (energies_.size() != size())
    throw std::runtime_error(
        "The number of energies does not match the number of structures in this molecular trajectory.");
  structureVector_.push_back(std::move(p));
  energies_.push_back(e);
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
  for (auto& s : structureVector_)
    s *= f;
  return *this;
}

const MolecularTrajectory& MolecularTrajectory::operator/=(double f) {
  for (auto& s : structureVector_)
    s /= f;
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
  bool hasSameSizeAsPreviousElementTypeCollection = ec.size() == molecularSize();
  bool noPositionsArePresent = empty();
  bool hasSameSizeAsPresentPositions = empty() || (structureVector_.front().size() == ec.size());
  return hasSameSizeAsPreviousElementTypeCollection || noPositionsArePresent || hasSameSizeAsPresentPositions;
}

bool MolecularTrajectory::additionOfPositionCollectionIsAllowed(const PositionCollection& p) const {
  bool validWithRespectToOtherPositionCollections = empty() || p.size() == structureVector_.front().size();
  bool validWithRespectToElementTypeCollection = p.rows() == molecularSize();

  bool validGivenThatElementCollectionIsAlreadySet = validWithRespectToElementTypeCollection;
  bool validGivenThatElementCollectionIsNotSet = (molecularSize() == 0 && validWithRespectToOtherPositionCollections);

  return validGivenThatElementCollectionIsAlreadySet || validGivenThatElementCollectionIsNotSet;
}

} /* namespace Utils */
} /* namespace Scine */
