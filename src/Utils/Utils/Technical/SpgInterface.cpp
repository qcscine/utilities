/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "SpgInterface.h"
#include "Utils/Geometry.h"
#include <spglib.h>
#include <map>
#include <unordered_set>
#define MAX_SYMMETRY_OPERATIONS 256

namespace Scine {
namespace Utils {

SpgInterface::CppCell SpgInterface::systemToCppCell(const PeriodicSystem& system, bool solidStateOnly) {
  const int nAtoms = system.atoms.size();
  int nResult = (solidStateOnly) ? static_cast<int>(system.solidStateAtomIndices.size()) : nAtoms;
  PositionCollection positions = PositionCollection::Zero(nResult, 3);
  auto endIter = system.solidStateAtomIndices.end();
  int shift = 0;
  std::vector<int> types;
  for (int i = 0; i < nAtoms; ++i) {
    if (solidStateOnly && system.solidStateAtomIndices.find(i) == endIter) {
      shift++;
      continue;
    }
    int index = i - shift;
    types.push_back(static_cast<int>(Utils::ElementInfo::Z(system.atoms[i].getElementType())));
    positions.row(index) = system.atoms.getPosition(i);
  }
  return {system.pbc, positions, types};
}

SpgInterface::CCell SpgInterface::systemToCell(const PeriodicSystem& system, bool solidStateOnly) {
  auto cppCell = systemToCppCell(system, solidStateOnly);
  return cppCellToCell(cppCell);
}

SpgInterface::CppCell SpgInterface::cellToCppCell(const SpgInterface::CCell& cell) {
  const auto [lattice, positions, types, nAtoms] = cell;
  // double array to Pbc
  Eigen::Matrix3d matrix = Eigen::Matrix3d::Zero();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      matrix(i, j) = lattice[i][j];
    }
  }
  matrix.transposeInPlace();
  const auto cellPbc = PeriodicBoundaries(matrix);

  // double array to PositionCollection
  const auto* cpos = *positions;
  PositionCollection relPositions = PositionCollection::Zero(nAtoms, 3);
  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < 3; ++j) {
      relPositions(i, j) = cpos[i][j];
    }
  }
  const PositionCollection absolutePositions = cellPbc.transform(relPositions);
  // array to vector
  std::vector<int> returnTypes;
  returnTypes.assign(*types, *types + nAtoms);
  return {cellPbc, absolutePositions, returnTypes};
}

SpgInterface::CCell SpgInterface::cppCellToCell(const SpgInterface::CppCell& cell) {
  auto [pbc, absPositions, types] = cell;
  const Eigen::Matrix3d& m = pbc.getCellMatrix().transpose();
  double lattice[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      lattice[i][j] = m(i, j);
    }
  }
  PositionCollection relPos = pbc.transform(absPositions, false);
  int nAtoms = types.size();
  auto* cPositions = new double[nAtoms][3];
  auto* cTypes = new int[nAtoms];
  for (int i = 0; i < nAtoms; ++i) {
    cTypes[i] = types[i];
    for (int j = 0; j < 3; ++j) {
      cPositions[i][j] = relPos(i, j);
    }
  }
  return {lattice, std::make_shared<double(*)[3]>(cPositions), std::make_shared<int*>(cTypes), nAtoms};
}

SpgInterface::CppCell SpgInterface::findPrimitiveCell(const PeriodicSystem& system, double epsilon, bool solidStateOnly) {
  auto cell = systemToCell(system, solidStateOnly);
  auto [lattice, positions, types, nAtoms] = cell;
  int nPrimitiveAtoms = spg_standardize_cell(lattice, *positions, *types, nAtoms, 1, 1, epsilon);
  if (nPrimitiveAtoms == 0) {
    auto error = spg_get_error_code();
    auto errorMsg = std::string(spg_get_error_message(error));
    throw std::runtime_error("Primitive cell could not be deduced\n" + errorMsg);
  }
  SpgInterface::CCell newCell(lattice, positions, types, nPrimitiveAtoms);
  return cellToCppCell(newCell);
}

PeriodicSystem SpgInterface::findPrimitiveCellSystem(const PeriodicSystem& system, double epsilon, bool solidStateOnly) {
  const auto cell = findPrimitiveCell(system, epsilon, solidStateOnly);
  const auto [primPbc, primPositions, primTypes] = cell;
  const auto primSize = primTypes.size();
  std::vector<ElementType> primElements;
  for (const auto& z : primTypes) {
    primElements.push_back(ElementType(z));
  }
  if (solidStateOnly) {
    // gather all indices
    std::unordered_set<unsigned> indices;
    for (unsigned i = 0; i < primSize; ++i) {
      indices.insert(i);
    }
    return PeriodicSystem(primPbc, primElements, primPositions, indices);
  }
  if (system.solidStateAtomIndices.empty()) {
    // we have no solid state indices
    return PeriodicSystem(primPbc, primElements, primPositions);
  }
  // have to determine new solid state indices
  // we first fold the old whole system into the primitive cell
  auto foldedSystem = PeriodicSystem(primPbc, system.atoms, system.solidStateAtomIndices);
  // we want check over all atoms of the folded system and check if we find the equivalent in the primitive system
  // this allows to map the indices
  // however SPGLib bugs out if positions are overlapping
  // get rid of them, while checking if they are of the same category + householding indices
  std::unordered_set<int> indicesToDiscard;
  const int nAtoms = foldedSystem.atoms.size();
  const double epsSquared = epsilon * epsilon;
  const std::vector<int> empty;
  for (int i = 0; i < nAtoms - 1; ++i) {
    auto iFind = foldedSystem.solidStateAtomIndices.find(i) == foldedSystem.solidStateAtomIndices.end();
    for (int j = i + 1; j < nAtoms; ++j) {
      if (Geometry::Distances::distanceSquared(foldedSystem.atoms.getPosition(i), foldedSystem.atoms.getPosition(j),
                                               foldedSystem.pbc) < epsSquared) {
        if (iFind != (foldedSystem.solidStateAtomIndices.find(j) == foldedSystem.solidStateAtomIndices.end())) {
          // makes sure we are not loosing solid state info
          throw std::runtime_error("Solid state and non-solid state atoms occupy the same position in the"
                                   "folded primitive cell");
        }
        indicesToDiscard.insert(j);
      }
    }
  }
  // remove atoms from foldedSystem
  int nNonClashing = nAtoms - indicesToDiscard.size();
  PositionCollection nonClashingPositions = PositionCollection(nNonClashing, 3);
  ElementTypeCollection nonClashingElements;
  int newIndex = 0;
  std::vector<int> discardedIndexToOriginalMap;
  for (int i = 0; i < nAtoms; ++i) {
    if (indicesToDiscard.find(i) == indicesToDiscard.end()) {
      nonClashingElements.push_back(foldedSystem.atoms.getElement(i));
      nonClashingPositions.row(newIndex) = foldedSystem.atoms.getPosition(i);
      discardedIndexToOriginalMap.push_back(i);
      newIndex++;
    }
  }
  foldedSystem.atoms.clear();
  foldedSystem.atoms.resize(nNonClashing);
  foldedSystem.atoms.setElements(nonClashingElements);
  foldedSystem.atoms.setPositions(nonClashingPositions);
  // use SPGLib symmetry operations so we don't miss any equalities
  auto foldedCell = systemToCppCell(foldedSystem);
  auto indexMap = mapSymmetryEquivalentPositions(cell, foldedCell, epsilon);
  std::unordered_set<unsigned> newIndices;
  auto endIter = system.solidStateAtomIndices.end();
  for (const auto& [primIndex, foldedIndices] : indexMap) {
    if (static_cast<int>(primSize) < primIndex) {
      throw std::runtime_error("Something went wrong while mapping the solid state atom indices");
    }
    if (foldedIndices.empty()) {
      throw std::runtime_error("Could not map every atom to an atom in the primitive cell");
    }
    for (const auto& index : foldedIndices) {
      int originalIndex = discardedIndexToOriginalMap[index];
      if (system.solidStateAtomIndices.find(originalIndex) != endIter) {
        // we found a solid state atoms
        newIndices.insert(primIndex);
        break;
      }
    }
  }
  return PeriodicSystem(primPbc, primElements, primPositions, newIndices);
}

std::map<int, std::vector<int>> SpgInterface::mapSymmetryEquivalentPositions(const CppCell& lhs, const CppCell& rhs,
                                                                             double epsilon) {
  std::map<int, std::vector<int>> result;
  const std::vector<int> empty;
  const auto lhsPositions = SpgInterface::getSymmetryEquivalents(lhs, epsilon);
  const auto rhsPositions = SpgInterface::getSymmetryEquivalents(rhs, epsilon);

  const int nAtomsLhs = static_cast<int>(lhs.types.size());
  const double epsSquared = epsilon * epsilon;
  for (const auto& pos : lhsPositions) {
    for (const auto& posRhs : rhsPositions) {
      // loop over lhs atoms
      for (int i = 0; i < nAtomsLhs; ++i) {
        auto [minDistance, minIndex] = minDistanceAndIndex(lhs.types[i], pos.row(i), rhs.types, posRhs, lhs.lattice);
        if (minDistance < epsSquared) {
          if (result.find(i) == result.end()) {
            result.insert(std::make_pair(i, empty));
          }
          result.at(i).push_back(minIndex);
        }
      }
    }
  }
  return result;
}

std::pair<int, int> SpgInterface::minDistanceAndIndex(int lhsType, const Position& lhsPos, std::vector<int> rhsTypes,
                                                      const PositionCollection& rhsPositions, const PeriodicBoundaries& pbc) {
  double minDistance = std::numeric_limits<double>::max();
  int minIndex = -1;
  int nAtomsRhs = static_cast<int>(rhsPositions.rows());
  // loop over rhs atoms
  for (int j = 0; j < nAtomsRhs; ++j) {
    if (lhsType != rhsTypes[j]) {
      continue;
    }
    const double dist = Utils::Geometry::Distances::distanceSquared(lhsPos, rhsPositions.row(j), pbc);
    if (dist < minDistance) {
      minDistance = dist;
      minIndex = j;
    }
  }
  return std::make_pair(minDistance, minIndex);
}

std::vector<SpgInterface::SymmetryOperation> SpgInterface::findSymmetryOperations(const CppCell& cell, double epsilon) {
  return findSymmetryOperations(cppCellToCell(cell), epsilon, cell.lattice.getCellMatrix());
}

std::vector<SpgInterface::SymmetryOperation> SpgInterface::findSymmetryOperations(const CCell& cell, double epsilon) {
  const auto& m = cell.lattice;
  Eigen::Matrix3d cellMatrix;
  // clang-format off
  cellMatrix << m[0][0], m[0][1], m[0][2],
                m[1][0], m[1][1], m[1][2],
                m[2][0], m[2][1], m[2][2];
  // clang-format on
  return findSymmetryOperations(cell, epsilon, cellMatrix.transpose());
}

std::vector<SpgInterface::SymmetryOperation> SpgInterface::findSymmetryOperations(const CCell& cell, double epsilon,
                                                                                  const Eigen::Matrix3d& pbc) {
  auto [lattice, positions, types, nAtoms] = cell;
  int rotations[MAX_SYMMETRY_OPERATIONS][3][3];
  double translations[MAX_SYMMETRY_OPERATIONS][3];
  int n_results =
      spg_get_symmetry(rotations, translations, MAX_SYMMETRY_OPERATIONS, lattice, *positions, *types, nAtoms, epsilon);
  if (n_results == 0) {
    auto error = spg_get_error_code();
    auto errorMsg = std::string(spg_get_error_message(error));
    throw std::runtime_error("Symmetry could not be deduced\n" + errorMsg);
  }
  std::vector<SymmetryOperation> result;
  result.reserve(n_results);
  for (int i = 0; i < n_results; ++i) {
    result.emplace_back(rotations[i], translations[i], pbc);
  }
  return result;
}

std::vector<PositionCollection> SpgInterface::getSymmetryEquivalents(const PositionCollection& positions,
                                                                     const std::vector<SpgInterface::SymmetryOperation>& operations) {
  std::vector<PositionCollection> result;
  result.reserve(operations.size());
  for (const auto& operation : operations) {
    PositionCollection pos = positions;
    applySymmetryOperation(pos, operation);
    result.push_back(pos);
  }
  return result;
}

std::vector<PositionCollection> SpgInterface::getSymmetryEquivalents(const CppCell& cell, double epsilon) {
  auto operations = findSymmetryOperations(cell, epsilon);
  return getSymmetryEquivalents(cell.positions, operations);
}

void SpgInterface::applySymmetryOperation(PositionCollection& positions, const SpgInterface::SymmetryOperation& operation) {
  positions *= operation.pbc.inverse(); // to relative coordinates
  using TransposedPC = Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor>;
  TransposedPC tmp1 = positions.transpose();
  TransposedPC tmp2 = operation.rotation * tmp1;
  positions = tmp2.transpose();
  Geometry::Manipulations::translatePositionsInPlace(positions, operation.translation);
  positions *= operation.pbc; // back to absolute Cartesian positions
}

std::vector<Eigen::Matrix3d> SpgInterface::findAlternativePbcs(const Eigen::Matrix3d& cellMatrix, double epsilon) {
  // we must NOT rely on CppCell / PeriodicBoundaries here, because this is meant to be used by the Pbc constructor
  double lattice[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      lattice[i][j] = cellMatrix(j, i); // shifted indices are on purpose
    }
  }
  // dummy position
  auto* cPositions = new double[1][3];
  cPositions[0][0] = 0.0;
  cPositions[0][1] = 0.0;
  cPositions[0][2] = 0.0;
  int cTypes[1] = {0};
  auto cell = CCell(lattice, std::make_shared<double(*)[3]>(cPositions), std::make_shared<int*>(cTypes), 1);
  auto operations = findSymmetryOperations(cell, epsilon);
  // translate cell matrix to PositionCollection to apply our method
  PositionCollection posOverload = Eigen::Map<const PositionCollection>(cellMatrix.data(), 3, 3);
  auto posResult = getSymmetryEquivalents(posOverload, operations);
  // translate back
  std::vector<Eigen::Matrix3d> result;
  result.reserve(posResult.size());
  for (const auto& pos : posResult) {
    Eigen::Matrix3d r = Eigen::Map<const Eigen::Matrix3d>(pos.data(), 3, 3);
    result.push_back(r);
  }
  return result;
}

bool SpgInterface::CppCell::isApprox(SpgInterface::CppCell rhs, double epsilon) const {
  // have to have the same lattice
  if (!lattice.isApprox(rhs.lattice, epsilon)) {
    return false;
  }
  // fast comparison in case they are equivalent without symmetry considerations
  if (positions.isApprox(rhs.positions, epsilon) && types == rhs.types) {
    return true;
  }
  const auto n = static_cast<int>(types.size());
  // check for constant translation
  Eigen::RowVector3d shift = Eigen::RowVector3d::Zero();
  for (int i = 0; i < n; ++i) {
    auto [dist, index] = SpgInterface::minDistanceAndIndex(types[i], positions.row(i), rhs.types, rhs.positions, lattice);
    if (dist > epsilon) {
      shift = positions.row(i) - rhs.positions.row(index);
      break;
    }
  }
  Geometry::Manipulations::translatePositionsInPlace(rhs.positions, shift);
  // fast comparison in case they were only constantly shifted
  if (positions.isApprox(rhs.positions, epsilon) && types == rhs.types) {
    return true;
  }
  const auto lhsPositions = SpgInterface::getSymmetryEquivalents(*this, epsilon);
  const auto rhsPositions = SpgInterface::getSymmetryEquivalents(rhs, epsilon);
  return isApproxImpl(rhs, epsilon, lhsPositions, rhsPositions);
}

bool SpgInterface::CppCell::isApprox(SpgInterface::CppCell rhs, double epsilon,
                                     const std::vector<SymmetryOperation>& lhsOperations,
                                     const std::vector<SymmetryOperation>& rhsOperations) const {
  // have to have the same lattice
  if (!lattice.isApprox(rhs.lattice, epsilon)) {
    return false;
  }
  // fast comparison in case they are equivalent without symmetry considerations
  if (positions.isApprox(rhs.positions, epsilon) && types == rhs.types) {
    return true;
  }
  const auto n = static_cast<int>(types.size());
  // check for constant translation
  Eigen::RowVector3d shift = Eigen::RowVector3d::Zero();
  for (int i = 0; i < n; ++i) {
    auto [dist, index] = SpgInterface::minDistanceAndIndex(types[i], positions.row(i), rhs.types, rhs.positions, lattice);
    if (dist > epsilon) {
      shift = positions.row(i) - rhs.positions.row(index);
      break;
    }
  }
  Geometry::Manipulations::translatePositionsInPlace(rhs.positions, shift);
  // fast comparison in case they were only constantly shifted
  if (positions.isApprox(rhs.positions, epsilon) && types == rhs.types) {
    return true;
  }
  const auto lhsPositions = SpgInterface::getSymmetryEquivalents(this->positions, lhsOperations);
  const auto rhsPositions = SpgInterface::getSymmetryEquivalents(rhs.positions, rhsOperations);
  return isApproxImpl(rhs, epsilon, lhsPositions, rhsPositions);
}

bool SpgInterface::CppCell::isApproxImpl(const SpgInterface::CppCell& rhs, double epsilon,
                                         const std::vector<PositionCollection>& lhsPositions,
                                         const std::vector<PositionCollection>& rhsPositions) const {
  const double epsSquared = epsilon * epsilon;
  const int nAtomsLhs = static_cast<int>(this->types.size());
  bool foundMatchingPositions = true;
  for (const auto& pos : lhsPositions) {
    for (const auto& posRhs : rhsPositions) {
      // loop over self
      foundMatchingPositions = true;
      for (int i = 0; i < nAtomsLhs; ++i) {
        auto [minDistance, _] = minDistanceAndIndex(types[i], pos.row(i), rhs.types, posRhs, lattice);
        if (minDistance > epsSquared) {
          foundMatchingPositions = false;
          break;
        }
      }
      if (foundMatchingPositions)
        break;
    }
    if (foundMatchingPositions)
      break;
  }
  return foundMatchingPositions;
}

} /* namespace Utils */
} /* namespace Scine */
