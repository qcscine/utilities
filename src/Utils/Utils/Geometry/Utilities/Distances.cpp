/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Distances.h"
#include "Properties.h"
#include <iomanip>

namespace Scine {
namespace Utils {
namespace Geometry {
namespace Distances {

double distanceSquared(const PositionCollection& p1, const PositionCollection& p2) {
  if (p1.rows() != p2.rows())
    throw std::runtime_error("The given PositionCollections have different number of positions.");
  double result = 0.0;
  for (int i = 0; i < p1.rows(); ++i) {
    Position pos1 = p1.row(i);
    Position pos2 = p2.row(i);
    result += distanceSquared(pos1, pos2);
  }
  return result;
}

double distanceSquared(const PositionCollection& p1, const PositionCollection& p2, const PeriodicBoundaries& pbc) {
  if (p1.rows() != p2.rows())
    throw std::runtime_error("The given PositionCollections have different number of positions.");
  double result = 0.0;
  for (int i = 0; i < p1.rows(); ++i) {
    Position pos1 = p1.row(i);
    Position pos2 = p2.row(i);
    result += distanceSquared(pos1, pos2, pbc);
  }
  return result;
}

double distanceSquared(const Position& p1, const Position& p2, const PeriodicBoundaries& pbc) {
  Position pos1 = pbc.translatePositionsIntoCell(p1);
  Position pos2 = pbc.translatePositionsIntoCell(p2);
  /* If real space distance smaller than half of smallest perpendicular, we can use fast algorithm */
  /* See "The Minimum Image Convention in Non-Cubic MD Cells" by W. Smith, 1989 for reference */
  if ((pos1 - pos2).squaredNorm() < 0.5 * pbc.getSmallestPerpendicularSquared())
    return pbc.fastMinimumImageDistanceSquared(pos1, pos2);
  else
    return pbc.bruteForceMinimumImageDistanceSquared(pos1, pos2);
}

BondOrderCollection nearestNeighborsBondOrders(const PositionCollection& positions, double margin) {
  BondOrderCollection bo = BondOrderCollection(positions.rows());
  for (int i = 0; i < positions.rows(); ++i) {
    Position pos = positions.row(i);
    std::vector<int> neighbors = nearestNeighborsInPositions(positions, pos, margin);
    for (const auto& neighbor : neighbors) {
      bo.setOrder(i, neighbor, 1.0);
    }
  }
  return bo;
}

BondOrderCollection nearestNeighborsBondOrders(const PositionCollection& positions, const PeriodicBoundaries& pbc,
                                               double margin) {
  BondOrderCollection bo = BondOrderCollection(positions.rows());
  for (int i = 0; i < positions.rows(); ++i) {
    Position pos = positions.row(i);
    std::vector<int> neighbors = nearestNeighborsInPositions(positions, pos, pbc, margin);
    for (const auto& neighbor : neighbors) {
      bo.setOrder(i, neighbor, 1.0);
    }
  }
  return bo;
}

std::vector<int> nearestNeighborsInPositions(const PositionCollection& positions, int index, double margin) {
  if (index >= positions.rows())
    throw std::runtime_error("The given index is out of range for the given PositionCollection.");
  Position pos = positions.row(index);
  return nearestNeighborsInPositions(positions, pos, margin);
}

std::vector<int> nearestNeighborsInPositions(const PositionCollection& positions, int index,
                                             const PeriodicBoundaries& pbc, double margin) {
  if (index >= positions.rows())
    throw std::runtime_error("The given index is out of range for the given PositionCollection.");
  Position pos = positions.row(index);
  return nearestNeighborsInPositions(positions, pos, pbc, margin);
}

std::vector<int> nearestNeighborsInPositions(const PositionCollection& positions, const Position& pos, double margin,
                                             double thresholdForSame) {
  std::map<int, double> preLimNeighborsMap;
  // init minimum distance with largest possible value but leave room for margin
  double minDistance = std::numeric_limits<double>::max() - margin - 1e-6;
  for (int i = 0; i < positions.rows(); ++i) {
    Position other = positions.row(i);
    double current = distance(pos, other);
    if (current < thresholdForSame) // the position in the Collection is identical with the given one -> ignore
      continue;
    if (current < minDistance + margin) {   // the position is a nearest neighbor
      if (current > minDistance - margin) { // it is in the same shell as the already known one
        preLimNeighborsMap.insert(std::make_pair(i, current));
      }
      else { // it must be closer than previous neighbor(s) beyond margin -> empty list and insert current one
        preLimNeighborsMap.clear();
        preLimNeighborsMap.insert(std::make_pair(i, current));
      }
    }
    if (current < minDistance) { // always want lowest possible minDistance of current shell
      minDistance = current;
    }
  }
  // Now check if really all neighbors in map qualify to be within margin
  // since they could have been added when the minimum distance was slightly longer
  std::vector<int> result;
  for (const auto& neighbor : preLimNeighborsMap) {
    if (neighbor.second < minDistance + margin) {
      result.push_back(neighbor.first);
    }
  }
  return result;
}

std::vector<int> nearestNeighborsInPositions(const PositionCollection& positions, const Position& pos,
                                             const PeriodicBoundaries& pbc, double margin, double thresholdForSame) {
  std::map<int, double> preLimNeighborsMap;
  // init minimum distance with largest possible value but leave room for margin
  double minDistance = std::numeric_limits<double>::max() - margin - 1e-6;
  for (int i = 0; i < positions.rows(); ++i) {
    Position other = positions.row(i);
    double current = distance(pos, other, pbc);
    if (current < thresholdForSame) // the position in the Collection is identical with the given one -> ignore
      continue;
    if (current < minDistance + margin) {   // the position is a nearest neighbor
      if (current > minDistance - margin) { // it is in the same shell as the already known one
        preLimNeighborsMap.insert(std::make_pair(i, current));
      }
      else { // it must be closer than previous neighbor(s) beyond margin -> empty list and insert current one
        preLimNeighborsMap.clear();
        preLimNeighborsMap.insert(std::make_pair(i, current));
      }
    }
    if (current < minDistance) { // always want lowest possible minDistance of current shell
      minDistance = current;
    }
  }
  // Now check if really all neighbors in map qualify to be within margin
  // since they could have been added when the minimum distance was slightly longer
  std::vector<int> result;
  for (const auto& neighbor : preLimNeighborsMap) {
    if (neighbor.second < minDistance + margin) {
      result.push_back(neighbor.first);
    }
  }
  return result;
}

std::vector<int> countAllNearestNeighbors(const PositionCollection& positions, double margin) {
  std::vector<int> result;
  result.reserve(positions.rows());
  for (int i = 0; i < positions.rows(); ++i) {
    result.push_back(countNearestNeighbors(positions, i, margin));
  }
  return result;
}

std::vector<int> countAllNearestNeighbors(const PositionCollection& positions, const PeriodicBoundaries& pbc, double margin) {
  std::vector<int> result;
  result.reserve(positions.rows());
  for (int i = 0; i < positions.rows(); ++i) {
    result.push_back(countNearestNeighbors(positions, i, pbc, margin));
  }
  return result;
}

int countNearestNeighbors(const PositionCollection& positions, int index, double margin) {
  if (index >= positions.rows())
    throw std::runtime_error("The given index is out of range for the given PositionCollection.");
  Position pos = positions.row(index);
  return countNearestNeighbors(positions, pos, margin);
}

int countNearestNeighbors(const PositionCollection& positions, int index, const PeriodicBoundaries& pbc, double margin) {
  if (index >= positions.rows())
    throw std::runtime_error("The given index is out of range for the given PositionCollection.");
  Position pos = positions.row(index);
  return countNearestNeighbors(positions, pos, pbc, margin);
}

int countNearestNeighbors(const PositionCollection& positions, const Position& pos, double margin, double thresholdForSame) {
  return nearestNeighborsInPositions(positions, pos, margin, thresholdForSame).size();
}

int countNearestNeighbors(PositionCollection positions, Position pos, const PeriodicBoundaries& pbc, double margin,
                          double thresholdForSame) {
  // more complex with PBC, because images could be within same shell
  // algorithm very similar to determine if nearest neighbor, but consider every image of every other position
  std::vector<double> preLimNeighborDistances;
  // init minimum distance with largest possible value but leave room for margin
  double minDistance = std::numeric_limits<double>::max() - margin - 1e-6;
  // we need to ensure that positions are within cell for proper allImageDistances
  pbc.translatePositionsIntoCellInPlace(positions);
  pbc.translatePositionsIntoCellInPlace(pos);
  for (int i = 0; i < positions.rows(); ++i) {
    Position other = positions.row(i);
    std::vector<double> distances = pbc.getAllImageDistancesSquared(pos, other);
    for (const auto& distance : distances) {
      double d = std::sqrt(distance);
      if (d < thresholdForSame)
        continue;                       // the position in the Collection is identical with the given one -> ignore
      if (d < minDistance + margin) {   // the position is a nearest neighbor
        if (d > minDistance - margin) { // it is in the same shell as the already known one
          preLimNeighborDistances.push_back(d);
        }
        else { // it must be closer than previous neighbor(s) beyond margin -> empty list and insert current one
          preLimNeighborDistances.clear();
          preLimNeighborDistances.push_back(d);
        }
      }
      if (d < minDistance) { // always want lowest possible minDistance of current shell
        minDistance = d;
      }
    }
  }
  // Now check if really all neighbors in map qualify to be within margin
  // since they could have been added when the minimum distance was slightly longer
  int counter = 0;
  for (const auto& distance : preLimNeighborDistances) {
    if (distance < minDistance + margin)
      counter++;
  }
  return counter;
}

int getIndexOfClosestAtom(const PositionCollection& positions, const Position& targetPosition,
                          double squaredDistanceConsideredZero) {
  assert(positions.rows() != 0 && "Cannot determine closest atom if there are no atoms!");
  double minimalDistanceSquared = std::numeric_limits<double>::max();
  int closestAtom = 0;

  auto nAtoms = static_cast<int>(positions.rows());

  for (int i = 0; i < nAtoms; i++) {
    double currentSquaredDistance = (positions.row(i) - targetPosition).squaredNorm();
    if (currentSquaredDistance <= squaredDistanceConsideredZero) {
      continue;
    }
    if (currentSquaredDistance < minimalDistanceSquared) {
      minimalDistanceSquared = currentSquaredDistance;
      closestAtom = i;
    }
  }

  return closestAtom;
}

int getIndexOfClosestAtom(const PositionCollection& positions, const Position& targetPosition,
                          const PeriodicBoundaries& pbc, double squaredDistanceConsideredZero) {
  assert(positions.rows() != 0 && "Cannot determine closest atom if there are no atoms!");
  double minimalDistanceSquared = std::numeric_limits<double>::max();
  int closestAtom = 0;

  auto nAtoms = static_cast<int>(positions.rows());

  for (int i = 0; i < nAtoms; i++) {
    Position currentPosition = positions.row(i);
    double currentSquaredDistance = distanceSquared(currentPosition, targetPosition, pbc);
    if (currentSquaredDistance <= squaredDistanceConsideredZero) {
      continue;
    }
    if (currentSquaredDistance < minimalDistanceSquared) {
      minimalDistanceSquared = currentSquaredDistance;
      closestAtom = i;
    }
  }

  return closestAtom;
}

std::vector<int> getIndicesCloseToAtom(const PositionCollection& positions, int targetAtomIndex,
                                       double distanceThreshold, bool sameAtomIncluded, bool upperTriangle) {
  std::vector<int> closeAtomsSet;
  auto nAtoms = static_cast<int>(positions.rows());
  int initialIndex = upperTriangle ? (sameAtomIncluded ? targetAtomIndex : targetAtomIndex + 1) : 0;

  double minimalDistanceSquared = std::numeric_limits<double>::min();

  Position targetPosition = positions.row(targetAtomIndex);

  for (int i = initialIndex; i < nAtoms; ++i) {
    double currentSquaredDistance = (positions.row(i) - targetPosition).norm();
    if (currentSquaredDistance <= distanceThreshold) {
      if (!sameAtomIncluded && currentSquaredDistance <= minimalDistanceSquared)
        continue;
      closeAtomsSet.push_back(i);
    }
  }
  return closeAtomsSet;
}

std::map<int, std::vector<int>> constructAtomPairList(const PositionCollection& positions, double distanceThreshold,
                                                      bool sameAtomIncluded, bool upperTriangle) {
  auto nAtoms = static_cast<int>(positions.rows());
  std::map<int, std::vector<int>> atomPairList;

  for (int atom = 0; atom < nAtoms; ++atom) {
    atomPairList.insert({atom, getIndicesCloseToAtom(positions, atom, distanceThreshold, sameAtomIncluded, upperTriangle)});
  }
  return atomPairList;
}

int getIndexOfAtomInStructure(const AtomCollection& structure, const Atom& atom, double squaredDistanceConsideredZero) {
  int index = 0;
  ElementType element = atom.getElementType();
  for (const auto& a : structure) {
    if (a.getElementType() == element) { // First check whether the element types match
      double squaredDistance = (a.getPosition() - atom.getPosition()).squaredNorm();
      if (squaredDistance <= squaredDistanceConsideredZero) { // Check whether the atoms have basically the same position
        return index;
      }
    }
    index++;
  }
  throw std::runtime_error("The given atom was not found in the given structure.");
}

std::vector<int> getListOfDivergingAtoms(const PositionCollection& reference, PositionCollection positions,
                                         double threshold, const ElementTypeCollection& elements) {
  std::vector<int> indices;
  indices.reserve(positions.rows());
  elements.empty() ? Manipulations::alignPositions(reference, positions)
                   : Manipulations::alignPositions(reference, positions, elements);
  Eigen::VectorXd difference = (reference - positions).rowwise().norm();
  for (int i = 0; i < difference.size(); ++i) {
    if (difference[i] > threshold) {
      indices.push_back(i);
    }
  }
  return indices;
}

std::vector<int> getListOfDivergingAtomsRobust(const PositionCollection& reference, PositionCollection positions,
                                               double threshold, double stopCriterion, int maxIterations,
                                               const ElementTypeCollection& elements, Core::Log log) {
  std::vector<int> indices;
  indices.reserve(positions.rows());
  Eigen::VectorXd rmsdVector = Eigen::VectorXd::Zero(positions.rows());
  Eigen::VectorXd oldRmsdVector;
  Eigen::VectorXd weightsStart =
      elements.empty()
          ? Eigen::VectorXd(Eigen::VectorXd::Ones(reference.rows()))
          : Eigen::VectorXd(Eigen::Map<const Eigen::VectorXd>(Properties::getMasses(elements).data(), elements.size()));
  Eigen::VectorXd weights = weightsStart;

  assert(maxIterations > 0);
  // Iterate
  log.output << std::setw(20) << "Iteration" << std::setw(20) << "Min RMSD" << std::setw(20) << "Max RMSD"
             << std::setw(20) << "Number Aligned" << Core::Log::nl;
  for (int iter = 0; iter < maxIterations; ++iter) {
    // Clear state at start of iter
    indices.clear();
    indices.reserve(positions.rows());
    oldRmsdVector = rmsdVector;

    // Align
    Manipulations::alignPositions(reference, positions, weights, rmsdVector);

    // Update weight and fill in indices
    for (int i = 0; i < rmsdVector.size(); ++i) {
      weights(i) = std::min(1. / rmsdVector(i), 20.);
      if (rmsdVector(i) > threshold)
        indices.push_back(i);
    }

    log.output << std::setw(20) << iter << std::setw(20) << rmsdVector.minCoeff() << std::setw(20)
               << rmsdVector.maxCoeff() << std::setw(20) << indices.size() << Core::Log::nl;
    // Check convergence
    if ((rmsdVector - oldRmsdVector).norm() < stopCriterion) {
      break;
    }
  }
  return indices;
}

} /* namespace Distances */
} // namespace Geometry
} // namespace Utils
} // namespace Scine
