/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "PeriodicSystem.h"
#include "AtomCollection.h"
#include "Utils/Bonds/BondDetector.h"
#include "Utils/Bonds/SolidStateBondDetector.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Geometry/GeometryUtilities.h"
#include "Utils/Technical/SpgInterface.h"
#include "Utils/Typenames.h"
#include <algorithm>
#include <sstream>

namespace Scine {
namespace Utils {

PeriodicSystem::PeriodicSystem(const PeriodicSystem& other)
  : PeriodicSystem(other.pbc, other.atoms, other.solidStateAtomIndices) {
}

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
  indicesCheck();
  canonicalize();
}

void PeriodicSystem::indicesCheck() const {
  if (!std::all_of(this->solidStateAtomIndices.begin(), this->solidStateAtomIndices.end(),
                   [&](unsigned i) { return static_cast<int>(i) < this->atoms.size(); })) {
    std::stringstream error;
    error << "At least one of the given solid state indices:\n";
    auto begin = solidStateAtomIndices.begin();
    const auto end = solidStateAtomIndices.end();
    // do not check for empty container, otherwise we wouldn't be here
    error << "[" << *begin++;
    while (begin != end) {
      error << ", " << *begin++;
    }
    error << "]\n";
    error << "is not valid for the given AtomCollection of size " << atoms.size();
    throw std::logic_error(error.str());
  }
}

void PeriodicSystem::canonicalize() {
  centerAndTranslateAtomsIntoCell();
  Eigen::Matrix3d rotation = pbc.getCanonicalizationRotationMatrix();
  if (rotation != Eigen::Matrix3d::Identity()) {
    this->pbc.canonicalize();
    this->atoms.setPositions(this->atoms.getPositions() * rotation);
  }
  centerAndTranslateAtomsIntoCell();
}

AtomCollection PeriodicSystem::getAtomCollectionWithImages(bool useSolidStateVanDerWaalsBonds) {
  return AtomCollection(atoms + getImageAtoms(useSolidStateVanDerWaalsBonds));
}

const AtomCollection& PeriodicSystem::getImageAtoms(bool useSolidStateVanDerWaalsBonds) {
  if (!_imageAtoms || _lastConstructionUsedVdwBonds != useSolidStateVanDerWaalsBonds) {
    constructImageAtoms(useSolidStateVanDerWaalsBonds);
  }
  return *_imageAtoms;
}

const std::unordered_map<unsigned, unsigned>& PeriodicSystem::getImageAtomsMap(bool useSolidStateVanDerWaalsBonds) {
  if (!_imageAtoms || _lastConstructionUsedVdwBonds != useSolidStateVanDerWaalsBonds) {
    constructImageAtoms(useSolidStateVanDerWaalsBonds);
  }
  return _imageAtomMap;
}

const BondOrderCollection& PeriodicSystem::getBondOrderCollectionWithImages(bool useSolidStateVanDerWaalsBonds) {
  if (!_imageBondOrderCollection || _lastConstructionUsedVdwBonds != useSolidStateVanDerWaalsBonds) {
    constructBondOrdersForImages(useSolidStateVanDerWaalsBonds);
  }
  return *_imageBondOrderCollection;
}

PeriodicSystem::masmData PeriodicSystem::getDataForMolassemblerInterpretation(bool useSolidStateVanDerWaalsBonds) {
  centerAndTranslateAtomsIntoCell();
  auto bondOrders = constructBondOrders(true, useSolidStateVanDerWaalsBonds);
  // remove 2nd shell if we rely on VdW, therefore pass boolean
  return getDataForMolassemblerInterpretation(bondOrders, useSolidStateVanDerWaalsBonds);
}

PeriodicSystem::masmData PeriodicSystem::getDataForMolassemblerInterpretation(const BondOrderCollection& bo,
                                                                              bool removeSolidSecondShell) {
  if (!_imageAtoms || atomsHaveChangedSinceLastImageConstruction()) {
    constructImageAtoms(bo, removeSolidSecondShell);
  }
  if (!_imageBondOrderCollection || atomsHaveChangedSinceLastImageConstruction()) {
    constructBondOrdersForImages(bo);
  }
  return std::make_tuple(atoms + *_imageAtoms, *_imageBondOrderCollection, solidStateAtomIndices, _imageAtomMap);
}

void PeriodicSystem::constructImageAtoms(bool useSolidStateVanDerWaalsBonds) {
  centerAndTranslateAtomsIntoCell();
  auto bondOrders = constructBondOrders(true, useSolidStateVanDerWaalsBonds);
  _lastConstructionUsedVdwBonds = useSolidStateVanDerWaalsBonds;
  constructImageAtoms(bondOrders, useSolidStateVanDerWaalsBonds); // remove 2nd shell if we rely on VdW
}

void PeriodicSystem::constructImageAtoms(const BondOrderCollection& bondOrders, bool removeSolidSecondShell) {
  clearImageAtoms();
  _imageAtoms = std::make_shared<AtomCollection>(AtomCollection());
  if (atoms.size() != bondOrders.getSystemSize()) {
    throw std::runtime_error("The atoms and bond orders do not match.");
  }
  auto itEnd = solidStateAtomIndices.end();
  const int n = atoms.size();
  for (int i = 0; i < n; ++i) {
    // because information currently only required if we want to remove second shell, we couple this bool to the look-up
    bool iIsSolid = removeSolidSecondShell && solidStateAtomIndices.find(i) != itEnd;
    for (int j = 0; j < i; ++j) {
      bool checkForSecondShell = iIsSolid && solidStateAtomIndices.find(j) != itEnd;
      if (bondOrders.getOrder(i, j) < 0.0) {
        Position p1 = atoms.getPosition(i);
        Position p2 = atoms.getPosition(j);
        auto diff = pbc.bruteForceMinimumImageDisplacementVector(p1, p2); // vector pointing from p1 to p2
        addPotentialImage(i, p2 - diff, checkForSecondShell);
        addPotentialImage(j, p1 + diff, checkForSecondShell);
      }
    }
  }
  _lastImageConstructedAtoms = atoms;
}

void PeriodicSystem::addPotentialImage(int index, const Position& position, bool checkForSecondShell) {
  if (_imageAtoms->size() > 0) {
    // if already images present, check whether we already created an image at this position
    int closestIndex = Geometry::Distances::getIndexOfClosestAtom(_imageAtoms->getPositions(), position);
    if (Geometry::Distances::distanceSquared(_imageAtoms->getPosition(closestIndex), position) < 0.01) {
      return;
    }
    if (checkForSecondShell) {
      // we rely here on nearest neighbors without periodic boundaries
      // if all neighbors of this potential image atom are other image atoms, it is part of the second shell
      PositionCollection allPosSoFar = PositionCollection(atoms.size() + _imageAtoms->size(), 3);
      allPosSoFar << atoms.getPositions(), _imageAtoms->getPositions();
      std::vector<int> neighbors = Geometry::Distances::nearestNeighborsInPositions(allPosSoFar, position);
      auto itEnd = _imageAtomMap.end();
      if (std::all_of(neighbors.begin(), neighbors.end(),
                      [&](int n) { return _imageAtomMap.find(static_cast<unsigned>(n)) != itEnd; })) {
        // all neighbors are in the imageAtomMap, therefore image atoms, so we assume second shell and don't add
        return;
      }
    }
  }
  _imageAtomMap.insert(std::make_pair(static_cast<unsigned>(atoms.size() + _imageAtoms->size()), static_cast<unsigned>(index)));
  _imageAtoms->push_back(Atom(atoms.getElement(index), position));
}

BondOrderCollection PeriodicSystem::constructBondOrders(bool periodic, bool useSolidStateVanDerWaalsBonds) const {
  indicesCheck();
  if (solidStateAtomIndices.empty()) {
    if (periodic) {
      return BondDetector::detectBonds(atoms, pbc, true);
    }
    return BondDetector::detectBonds(atoms);
  }
  if (periodic) {
    return SolidStateBondDetector::detectBonds(atoms, pbc, solidStateAtomIndices, true, useSolidStateVanDerWaalsBonds);
  }
  return SolidStateBondDetector::detectBonds(atoms, solidStateAtomIndices, useSolidStateVanDerWaalsBonds);
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

void PeriodicSystem::constructBondOrdersForImages(bool useSolidStateVanDerWaalsBonds) {
  centerAndTranslateAtomsIntoCell();
  auto bondOrders = constructBondOrders(true, useSolidStateVanDerWaalsBonds);
  _lastConstructionUsedVdwBonds = useSolidStateVanDerWaalsBonds;
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

PeriodicSystem PeriodicSystem::getPrimitiveCellSystem(double epsilon, bool solidStateOnly) const {
  return SpgInterface::findPrimitiveCellSystem(*this, epsilon, solidStateOnly);
}

bool PeriodicSystem::isApprox(const PeriodicSystem& other, double eps) const {
  const auto nAtoms = atoms.size();
  if (nAtoms != other.atoms.size()) {
    return false;
  }
  const auto nIndices = solidStateAtomIndices.size();
  if (!pbc.isApprox(other.pbc, eps) || nIndices != other.solidStateAtomIndices.size()) {
    return false;
  }
  // fast position comparison
  if (atoms.isApprox(other.atoms, eps)) {
    return true;
  }
  // slow position comparison taking periodicity into account
  // naively compare primitive cell
  const auto prim = SpgInterface::findPrimitiveCell(*this, eps);
  const auto primRhs = SpgInterface::findPrimitiveCell(other, eps);
  if (prim.isApprox(primRhs, eps)) {
    return true;
  }
  // compare primitive cell solid state indices only
  auto primSolid = SpgInterface::findPrimitiveCell(*this, eps, true);
  auto primSolidRhs = SpgInterface::findPrimitiveCell(other, eps, true);
  if (!primSolid.isApprox(primSolidRhs, eps)) {
    return false;
  }

  // get symmetry operations
  auto lhsOperations = SpgInterface::findSymmetryOperations(primSolid, eps);
  auto rhsOperations = SpgInterface::findSymmetryOperations(primSolidRhs, eps);

  // get non-solid-state-atoms
  auto endIter = this->solidStateAtomIndices.end();
  auto endIterRhs = other.solidStateAtomIndices.end();
  PositionCollection positions = PositionCollection::Zero(static_cast<int>(nAtoms - nIndices), 3);
  PositionCollection positionsRhs = PositionCollection::Zero(static_cast<int>(nAtoms - nIndices), 3);
  std::vector<int> elements;
  std::vector<int> elementsRhs;
  int lhsIndex = 0;
  int rhsIndex = 0;
  for (int i = 0; i < nAtoms; ++i) {
    if (this->solidStateAtomIndices.find(i) == endIter) {
      positions.row(lhsIndex++) = this->atoms.getPosition(i);
      elements.push_back(static_cast<int>(ElementInfo::Z(this->atoms.getElement(i))));
    }
    if (other.solidStateAtomIndices.find(i) == endIterRhs) {
      positionsRhs.row(rhsIndex++) = other.atoms.getPosition(i);
      elementsRhs.push_back(static_cast<int>(ElementInfo::Z(other.atoms.getElement(i))));
    }
  }
  auto lhsCell = SpgInterface::CppCell(primSolid.lattice, positions, elements);
  auto rhsCell = SpgInterface::CppCell(primSolidRhs.lattice, positionsRhs, elementsRhs);
  return lhsCell.isApprox(rhsCell, eps, lhsOperations, rhsOperations);
}

bool PeriodicSystem::operator==(const PeriodicSystem& other) const {
  return isApprox(other, 1e-12);
}

bool PeriodicSystem::operator!=(const PeriodicSystem& other) const {
  return !(*this == other);
}

PeriodicSystem PeriodicSystem::operator*(const Eigen::Vector3i& scalingFactors) const {
  PeriodicSystem super = PeriodicSystem(*this);
  super *= scalingFactors;
  return super;
}

PeriodicSystem PeriodicSystem::operator*(const std::vector<int>& scalingFactors) const {
  PeriodicSystem super = PeriodicSystem(*this);
  super *= scalingFactors;
  return super;
}

PeriodicSystem& PeriodicSystem::operator*=(const std::vector<int>& scalingFactors) {
  if (scalingFactors.size() != 3) {
    throw std::runtime_error("Scaling factor of PeriodicSystem must be either a number or a vector of length 3");
  }
  const Eigen::Vector3i vec = Eigen::Map<const Eigen::Vector3i>(scalingFactors.data(), 3);
  this->operator*=(vec);
  return *this;
}

PeriodicSystem& PeriodicSystem::operator*=(const Eigen::Vector3i& scalingFactors) {
  if ((scalingFactors.array() < 1).any()) {
    throw std::runtime_error(
        "At least one scaling factor of the periodic system is less than 1, which is not possible.");
  }
  for (int i = 0; i < 3; ++i) {
    if (scalingFactors[i] > 1 && !pbc.getPeriodicity()[i]) {
      throw std::runtime_error(
          "Scaling factor of the periodic system is larger than 1 for a dimension which is not periodic.");
    }
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
  PeriodicSystem super = (*this);
  super *= scalingFactor;
  return super;
}

PeriodicSystem& PeriodicSystem::operator*=(int scalingFactor) {
  if (scalingFactor < 1) {
    throw std::runtime_error("Specified scaling factor of " + std::to_string(scalingFactor) + ", but it must be at least 1");
  }
  Eigen::Vector3i scaling = Eigen::Vector3i::Constant(scalingFactor);
  for (int i = 0; i < 3; ++i) {
    if (!pbc.getPeriodicity()[i]) {
      scaling[i] = 1;
    }
  }
  *this *= scaling;
  return *this;
}

} /* namespace Utils */
} /* namespace Scine */
