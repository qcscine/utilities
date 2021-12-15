/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "PeriodicSystem.h"
#include "Utils/Bonds/BondDetector.h"
#include "Utils/Bonds/SolidStateBondDetector.h"
#include "Utils/Geometry/GeometryUtilities.h"

namespace Scine {
namespace Utils {

PeriodicSystem::PeriodicSystem(const PeriodicBoundaries& pbc, int N, std::unordered_set<unsigned> solidStateAtomIndices)
  : PeriodicSystem(pbc, AtomCollection(N), std::move(solidStateAtomIndices)) {
}

PeriodicSystem::PeriodicSystem(const PeriodicBoundaries& pbc, const ElementTypeCollection& elements,
                               const PositionCollection& positions, std::unordered_set<unsigned> solidStateAtomIndices)
  : PeriodicSystem(pbc, AtomCollection(elements, positions), std::move(solidStateAtomIndices)) {
}

PeriodicSystem::PeriodicSystem(const PeriodicBoundaries& pbc, AtomCollection atoms,
                               std::unordered_set<unsigned> solidStateAtomIndices)
  : pbc(pbc), atoms(std::move(atoms)), solidStateAtomIndices(std::move(solidStateAtomIndices)) {
  if (!std::all_of(this->solidStateAtomIndices.begin(), this->solidStateAtomIndices.end(),
                   [&](unsigned i) { return static_cast<int>(i) < this->atoms.size(); })) {
    throw std::logic_error("At least one of the given solid state indices is not valid for the given AtomCollection.");
  }
  centerAndTranslateAtomsIntoCell();
}

AtomCollection PeriodicSystem::getAtomCollectionWithImages() {
  return AtomCollection(atoms + getImageAtoms());
}

const AtomCollection& PeriodicSystem::getImageAtoms() {
  if (!_imageAtoms) {
    constructImageAtoms();
  }
  return *_imageAtoms;
}

const std::unordered_map<unsigned, unsigned>& PeriodicSystem::getImageAtomsMap() {
  if (!_imageAtoms) {
    constructImageAtoms();
  }
  return _imageAtomMap;
}

const BondOrderCollection& PeriodicSystem::getBondOrderCollectionWithImages() {
  if (!_imageBondOrderCollection) {
    constructBondOrdersForImages();
  }
  return *_imageBondOrderCollection;
}

PeriodicSystem::masmData PeriodicSystem::getDataForMolassemblerInterpretation() {
  centerAndTranslateAtomsIntoCell();
  auto bondOrders = constructBondOrders();
  return getDataForMolassemblerInterpretation(bondOrders);
}

PeriodicSystem::masmData PeriodicSystem::getDataForMolassemblerInterpretation(const BondOrderCollection& bo) {
  if (!_imageAtoms || atomsHaveChangedSinceLastImageConstruction()) {
    constructImageAtoms(bo);
  }
  if (!_imageBondOrderCollection || atomsHaveChangedSinceLastImageConstruction()) {
    constructBondOrdersForImages(bo);
  }
  return std::make_tuple(atoms + *_imageAtoms, *_imageBondOrderCollection, solidStateAtomIndices, _imageAtomMap);
}

void PeriodicSystem::constructImageAtoms() {
  centerAndTranslateAtomsIntoCell();
  auto bondOrders = constructBondOrders();
  constructImageAtoms(bondOrders);
}

void PeriodicSystem::constructImageAtoms(const BondOrderCollection& bondOrders) {
  clearImageAtoms();
  _imageAtoms = std::make_shared<AtomCollection>(AtomCollection());
  if (atoms.size() != bondOrders.getSystemSize()) {
    throw std::runtime_error("The atoms and bond orders do not match.");
  }
  const int n = atoms.size();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      if (bondOrders.getOrder(i, j) < 0.0) {
        Position p1 = atoms.getPosition(i);
        Position p2 = atoms.getPosition(j);
        auto diff = pbc.bruteForceMinimumImageDisplacementVector(p1, p2); // vector pointing from p1 to p2
        addPotentialImage(i, p2 - diff);
        addPotentialImage(j, p1 + diff);
      }
    }
  }
  _lastImageConstructedAtoms = atoms;
}

void PeriodicSystem::addPotentialImage(int index, const Position& position) {
  if (_imageAtoms->size() > 0) {
    // if already images present, check whether we already created an image at this position
    int closestIndex = Geometry::Distances::getIndexOfClosestAtom(_imageAtoms->getPositions(), position);
    if (Geometry::Distances::distanceSquared(_imageAtoms->getPosition(closestIndex), position) < 0.01) {
      return;
    }
  }
  _imageAtomMap.insert(std::make_pair(static_cast<unsigned>(atoms.size() + _imageAtoms->size()), static_cast<unsigned>(index)));
  _imageAtoms->push_back(Atom(atoms.getElement(index), position));
}

BondOrderCollection PeriodicSystem::constructBondOrders(bool periodic) const {
  if (!std::all_of(solidStateAtomIndices.begin(), solidStateAtomIndices.end(), [&](int i) { return i < atoms.size(); })) {
    throw std::logic_error("At least one of the given solid state indices is not valid for the given AtomCollection.");
  }
  if (solidStateAtomIndices.empty()) {
    if (periodic) {
      return BondDetector::detectBonds(atoms, pbc, true);
    }
    return BondDetector::detectBonds(atoms);
  }
  if (periodic) {
    return SolidStateBondDetector::detectBonds(atoms, pbc, solidStateAtomIndices, true);
  }
  return SolidStateBondDetector::detectBonds(atoms, solidStateAtomIndices);
}

void PeriodicSystem::makeBondOrdersAcrossBoundariesNegative(BondOrderCollection& bondOrders) const {
  if (atoms.size() != bondOrders.getSystemSize()) {
    throw std::runtime_error("The atoms and bond orders do not match.");
  }
  bondOrders.setToAbsoluteValues();
  const int n = bondOrders.getSystemSize();
  for (int i = 0; i < n; ++i) {
    auto pos = atoms.getPosition(i);
    for (int j = 0; j < i; ++j) {
      double order = bondOrders.getOrder(i, j);
      if (order > 0.0 && pbc.minimumDistanceViaImage(pos, atoms.getPosition(j))) {
        bondOrders.setOrder(i, j, -1.0 * order);
      }
    }
  }
}

void PeriodicSystem::constructBondOrdersForImages() {
  centerAndTranslateAtomsIntoCell();
  auto bondOrders = constructBondOrders();
  constructBondOrdersForImages(bondOrders);
}

void PeriodicSystem::constructBondOrdersForImages(const BondOrderCollection& bondOrdersWithoutImages) {
  if (atoms.size() != bondOrdersWithoutImages.getSystemSize()) {
    throw std::runtime_error("The given bond orders without images do not fit the given atoms.");
  }
  if (!_imageAtoms || atomsHaveChangedSinceLastImageConstruction()) {
    constructImageAtoms(bondOrdersWithoutImages);
  }
  const unsigned n = atoms.size();
  const int N = atoms.size() + _imageAtoms->size();
  if (!_imageBondOrderCollection) {
    _imageBondOrderCollection = std::make_shared<BondOrderCollection>(BondOrderCollection());
  }
  _imageBondOrderCollection->resize(N);
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < i; ++j) {
      double order = bondOrdersWithoutImages.getOrder(i, j);
      _imageBondOrderCollection->setOrder(i, j, std::fabs(order));
      // if negative order set bonds with images based on map
      if (order < 0.0) {
        // search for value (real atom index) in map to make bond with key (image atom index)
        for (const auto& it : _imageAtomMap) {
          if (it.second == i) {
            _imageBondOrderCollection->setOrder(it.first, i, std::fabs(order));
          }
          if (it.second == j) {
            _imageBondOrderCollection->setOrder(it.first, j, std::fabs(order));
          }
        }
      }
    }
  }
}

void PeriodicSystem::centerAndTranslateAtomsIntoCell() {
  auto com = Geometry::Properties::getCenterOfMass(atoms);
  Position centerOfCell = Position::Constant(0.5);
  pbc.transformInPlace(centerOfCell);
  Displacement shift = centerOfCell - com;
  auto positions = atoms.getPositions();
  Geometry::Manipulations::translatePositionsInPlace(positions, shift);
  atoms.setPositions(positions);
  translateAtomsIntoCell();
}

void PeriodicSystem::translateAtomsIntoCell() {
  clearImageAtoms();
  auto positions = this->pbc.translatePositionsIntoCell(this->atoms.getPositions());
  this->atoms.setPositions(positions);
}

void PeriodicSystem::clear() {
  atoms.clear();
  clearImageAtoms();
}

bool PeriodicSystem::operator==(const PeriodicSystem& other) const {
  return (atoms == other.atoms && pbc == other.pbc && solidStateAtomIndices == other.solidStateAtomIndices);
}

bool PeriodicSystem::operator!=(const PeriodicSystem& other) const {
  return !(*this == other);
}

PeriodicSystem PeriodicSystem::operator*(const Eigen::Vector3i& scalingFactors) const {
  PeriodicSystem super = PeriodicSystem(this->pbc, this->atoms, this->solidStateAtomIndices);
  return super *= scalingFactors;
}

PeriodicSystem& PeriodicSystem::operator*=(const Eigen::Vector3i& scalingFactors) {
  if ((scalingFactors.array() < 1).any()) {
    throw std::runtime_error(
        "At least one scaling factor of the periodic system is less than 1, which is not possible.");
  }
  // do operations within cell -> save original state to translate back afterwards
  const PositionCollection originalPositions = atoms.getPositions();
  centerAndTranslateAtomsIntoCell();
  auto basisPositions = atoms.getPositions();
  const PositionCollection diff = basisPositions - originalPositions;
  const int nAtoms = atoms.size();
  // save basis data
  const std::unordered_set<unsigned> originalSolidStateIndices = solidStateAtomIndices;
  const ElementTypeCollection elements = atoms.getElements();
  pbc.transformInPlace(basisPositions, false);
  // go through all relative shifts and add another copy of the atoms to the atoms
  for (int a = 0; a < scalingFactors[0]; ++a) {
    for (int b = 0; b < scalingFactors[1]; ++b) {
      for (int c = 0; c < scalingFactors[2]; ++c) {
        Eigen::RowVector3d shiftVector;
        shiftVector << a, b, c;
        if (shiftVector.isApprox(Eigen::RowVector3d::Zero())) {
          continue; // skip 0 0 0, because we already have basis atoms
        }
        // add a to x, b to y and c to z
        PositionCollection newPositions = basisPositions.array().rowwise() + shiftVector.array();
        auto newAtoms = AtomCollection(elements, pbc.transform(newPositions));
        atoms += newAtoms;
        // add solid state indices for new atoms
        int indexShift = atoms.size() - nAtoms;
        for (auto index : originalSolidStateIndices) {
          solidStateAtomIndices.insert(index + indexShift);
        }
      }
    }
  }
  // shift positions back possibly outside of cell
  auto positions = atoms.getPositions();
  for (int i = 0; i < atoms.size(); i += nAtoms) {
    for (int j = 0; j < nAtoms; ++j) {
      positions.row(i + j) -= diff.row(j);
    }
  }
  atoms.setPositions(positions);
  // adjust periodic boundaries
  Eigen::Vector3d sf = scalingFactors.cast<double>();
  this->pbc *= sf;
  return *this;
}

PeriodicSystem PeriodicSystem::operator*(int scalingFactor) const {
  Eigen::Vector3i scaling = Eigen::Vector3i::Constant(scalingFactor);
  return *this * scaling;
}

PeriodicSystem& PeriodicSystem::operator*=(int scalingFactor) {
  Eigen::Vector3i scaling = Eigen::Vector3i::Constant(scalingFactor);
  return *this *= scaling;
}

} /* namespace Utils */
} /* namespace Scine */
