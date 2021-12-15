/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_GEOMETRY_DISTANCES_H_
#define UTILS_GEOMETRY_DISTANCES_H_

#include "Manipulations.h"
#include "Utils/DataStructures/PeriodicBoundaries.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Typenames.h"
#include <Core/Log.h>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace Scine {
namespace Utils {

/// @brief Functionalities working with an entire geometry (PositionCollection).
namespace Geometry {

/// @brief Functionalities calculating distances based on geometry or giving distance-based information.
namespace Distances {
/**
 * @brief Get the distance squared between two PositionCollections
 * @param p1 First PositionCollection
 * @param p2 Second PositionCollection
 * @return Squared distance
 */
double distanceSquared(const PositionCollection& p1, const PositionCollection& p2);

/**
 * @brief Get the distance squared between two PositionCollections using the minimum image convention
 * @param p1 First PositionCollection
 * @param p2 Second PositionCollection
 * @param pbc Periodic Boundary Condition
 * @return Squared distance
 */
double distanceSquared(const PositionCollection& p1, const PositionCollection& p2, const PeriodicBoundaries& pbc);

/**
 * @brief Get the distance squared between two Positions
 * @param p1 First Position
 * @param p2 Second Position
 * @return Squared distance
 */
inline double distanceSquared(const Position& p1, const Position& p2) {
  return (p1 - p2).squaredNorm();
}

/**
 * @brief Get the distance squared between two Positions using the minimum image convention
 * @param p1 First Position
 * @param p2 Second Position
 * @param pbc Periodic Boundary Condition
 * @return Squared distance
 */
double distanceSquared(const Position& p1, const Position& p2, const PeriodicBoundaries& pbc);

/**
 * @brief Get the distance between two PositionCollections
 * @param p1 First PositionCollection
 * @param p2 Second PositionCollection
 * @return Distance
 */
inline double distance(const PositionCollection& p1, const PositionCollection& p2) {
  return std::sqrt(distanceSquared(p1, p2));
}

/**
 * @brief Get the distance between two PositionCollections using the minimum image convention
 * @param p1 First PositionCollection
 * @param p2 Second PositionCollection
 * @param pbc Periodic Boundary Condition
 * @return Distance
 */
inline double distance(const PositionCollection& p1, const PositionCollection& p2, const PeriodicBoundaries& pbc) {
  return std::sqrt(distanceSquared(p1, p2, pbc));
}

/**
 * @brief Get the distance between two Positions
 * @param p1 First Position
 * @param p2 Second Position
 * @return Distance
 */
inline double distance(const Position& p1, const Position& p2) {
  return (p1 - p2).norm();
}

/**
 * @brief Get the distance between two Positions using the minimum image convention
 * @param p1 First Position
 * @param p2 Second Position
 * @param pbc Periodic Boundary Condition
 * @return Distance
 */
inline double distance(const Position& p1, const Position& p2, const PeriodicBoundaries& pbc) {
  return std::sqrt(distanceSquared(p1, p2, pbc));
}

/**
 * @brief Get the BondOrderCollection of a Geometry where only nearest neighbors are connected by single bonds.
 *
 * @note The criterion for a bond is that at least ONE of the two positions has to be a nearest neighbor of the other
 * one, and the nearest neighbor property does NOT have to be mutual. This means that the resulting bonding partner(s)
 * of a position can be different to its nearest neighbors by the position being a nearest neighbor of its bonding
 * partner(s).
 * in short: if x is NN of y OR y is a NN of x: bond between x and y
 * This also means that each position has at least one bond.
 *
 * @param positions The PositionCollection.
 * @param margin The margin to consider two distances to be equal
 * @return BondOrderCollection A matrix containing 1 and 0 for bonds with each nearest neighbor of each position.
 */
BondOrderCollection nearestNeighborsBondOrders(const PositionCollection& positions, double margin = 0.1);

/**
 * @brief Get the BondOrderCollection of a Geometry where only nearest neighbors are connected by single bonds with
 * periodic boundaries considered.
 *
 * @note The criterion for a bond is that at least ONE of the two positions has to be a nearest neighbor of the other
 * one, and the nearest neighbor property does NOT have to be mutual. This means that the resulting bonding partner(s)
 * of a position can be different to its nearest neighbors by the position being a nearest neighbor of its bonding
 * partner(s).
 * in short: if x is NN of y OR y is a NN of x: bond between x and y
 * This also means that each position has at least one bond.
 *
 * @param positions The PositionCollection.
 * @param pbc Periodic Boundary Condition.
 * @param margin The margin to consider two distances to be equal
 * @return BondOrderCollection A matrix containing 1 and 0 for bonds with each nearest neighbor of each position.
 */
BondOrderCollection nearestNeighborsBondOrders(const PositionCollection& positions, const PeriodicBoundaries& pbc,
                                               double margin = 0.1);

/**
 * @brief Get all indices of the nearest neighbors of one position with a PositionCollection.
 * @param positions The PositionCollection.
 * @param index The index of the Position you want the nearest neighbors from starting from 0.
 * @param margin The margin to consider two distances to be equal
 * @return std::vector<int> An array of indices that are nearest neighbors of the given position
 */
std::vector<int> nearestNeighborsInPositions(const PositionCollection& positions, int index, double margin = 0.1);

/**
 * @brief Get all indices of the nearest neighbors of one position with a PositionCollection with periodic boundaries
 * considered.
 * @param positions The PositionCollection.
 * @param index The index of the Position you want the nearest neighbors from starting from 0.
 * @param pbc Periodic Boundary Condition.
 * @param margin The margin to consider two distances to be equal
 * @return std::vector<int> An array of indices that are nearest neighbors of the given position
 */
std::vector<int> nearestNeighborsInPositions(const PositionCollection& positions, int index,
                                             const PeriodicBoundaries& pbc, double margin = 0.1);

/**
 * @brief Get all indices of the nearest neighbors of one position with a PositionCollection. The position can be within
 * the PositionCollection or not.
 * @param positions The PositionCollection.
 * @param pos The Position you want the nearest neighbors from.
 * @param margin The margin to consider two distances to be equal
 * @param thresholdForSame The distance under which the given Position is assumed to be identical to a Position within
 * the given collection and is therefore not considered a nearest neighbor.
 * @return std::vector<int> An array of indices that are nearest neighbors of the given position
 */
std::vector<int> nearestNeighborsInPositions(const PositionCollection& positions, const Position& pos,
                                             double margin = 0.1, double thresholdForSame = 0.01);

/**
 * @brief Get all indices of the nearest neighbors of one position with a PositionCollection with periodic boundaries
 * considered. The position can be within the PositionCollection or not.
 * @param positions The PositionCollection.
 * @param pos The Position you want the nearest neighbors from.
 * @param pbc Periodic Boundary Condition.
 * @param margin The margin to consider two distances to be equal
 * @param thresholdForSame The distance under which the given Position is assumed to be identical to a Position within
 * the given collection and is therefore not considered a nearest neighbor.
 * @return std::vector<int> An array of indices that are nearest neighbors of the given position
 */
std::vector<int> nearestNeighborsInPositions(const PositionCollection& positions, const Position& pos,
                                             const PeriodicBoundaries& pbc, double margin = 0.1,
                                             double thresholdForSame = 0.01);

/**
 * @brief Get the number of the nearest neighbors for all Positions within a PositionCollection.
 * @param positions The PositionCollection.
 * @param margin The margin to consider two distances to be equal
 * @return std::vector<int> An array with on entry per position giving the number of nearest neighbors this position has.
 */
std::vector<int> countAllNearestNeighbors(const PositionCollection& positions, double margin = 0.1);

/**
 * @brief Get the number of the nearest neighbors for all Positions within a PositionCollection with periodic boundaries
 * considered.
 * @param positions The PositionCollection.
 * @param pbc Periodic Boundary Condition.
 * @param margin The margin to consider two distances to be equal
 * @return std::vector<int> An array with on entry per position giving the number of nearest neighbors this position has.
 */
std::vector<int> countAllNearestNeighbors(const PositionCollection& positions, const PeriodicBoundaries& pbc,
                                          double margin = 0.1);

/**
 * @brief Get the number of nearest neighbors of one position with a PositionCollection.
 * @param positions The PositionCollection.
 * @param index The index of the Position you want the nearest neighbors from starting from 0.
 * @param margin The margin to consider two distances to be equal
 * @return int The number of nearest neighbors.
 */
int countNearestNeighbors(const PositionCollection& positions, int index, double margin = 0.1);

/**
 * @brief Get the number of nearest neighbors of one position with a PositionCollection with periodic boundaries
 * considered.
 * @param positions The PositionCollection.
 * @param index The index of the Position you want the nearest neighbors from starting from 0.
 * @param pbc Periodic Boundary Condition.
 * @param margin The margin to consider two distances to be equal
 * @return int The number of nearest neighbors.
 */
int countNearestNeighbors(const PositionCollection& positions, int index, const PeriodicBoundaries& pbc, double margin = 0.1);

/**
 * @brief Get the number of nearest neighbors of one position with a PositionCollection with periodic boundaries
 * considered. The position can be within the PositionCollection or not.
 * @param positions The PositionCollection.
 * @param pos The Position you want the nearest neighbors from.
 * @param pbc Periodic Boundary Condition.
 * @param margin The margin to consider two distances to be equal
 * @param thresholdForSame The distance under which the given Position is assumed to be identical to a Position within
 * the given collection and is therefore not considered a nearest neighbor.
 * @return int The number of nearest neighbors.
 */
int countNearestNeighbors(const PositionCollection& positions, const Position& pos, double margin = 0.1,
                          double thresholdForSame = 0.01);

/**
 * @brief Get the number of nearest neighbors of one position with a PositionCollection with periodic boundaries
 * considered. The position can be within the PositionCollection or not.
 * @param positions The PositionCollection.
 * @param pos The Position you want the nearest neighbors from.
 * @param pbc Periodic Boundary Condition.
 * @param margin The margin to consider two distances to be equal
 * @param thresholdForSame The distance under which the given Position is assumed to be identical to a Position within
 * the given collection and is therefore not considered a nearest neighbor.
 * @return int The number of nearest neighbors.
 */
int countNearestNeighbors(const PositionCollection& positions, const Position& pos, const PeriodicBoundaries& pbc,
                          double margin = 0.1, double thresholdForSame = 0.01);

/**
 * @brief Get the index of closest position (atom) to a given position in space.
 * @param positions A set of positions to be traversed.
 * @param targetPosition The target position.
 * @param distanceConsideredZero Squared distance between two positions resulting in them being considered equal
 *                               and therefore skipping this position as it is already present in the given
 *                               PositionCollection. The default is set to a negative value, so that all positions
 *                               are considered as possible candidates for the closest one.
 * @return int The index of the closest atom to the given target.
 */
int getIndexOfClosestAtom(const PositionCollection& positions, const Position& targetPosition,
                          double squaredDistanceConsideredZero = -1.0);

/**
 * @brief Get the index of closest position (atom) to a given position in space with periodic boundaries considered.
 * @param positions A set of positions to be traversed.
 * @param targetPosition The target position.
 * @param pbc The Periodic Boundary Conditions.
 * @param distanceConsideredZero Squared distance between two positions resulting in them being considered equal
 *                               and therefore skipping this position as it is already present in the given
 *                               PositionCollection. The default is set to a negative value, so that all positions
 *                               are considered as possible candidates for the closest one.
 * @return int The index of the closest atom to the given target.
 */
int getIndexOfClosestAtom(const PositionCollection& positions, const Position& targetPosition,
                          const PeriodicBoundaries& pbc, double squaredDistanceConsideredZero = -1.0);

/**
 * @brief Get the index of all the atoms within a certain distance to a given position in space.
 * @param positions A set of positions to be traversed.
 * @param targetAtomIndex The index of the target atom.
 * @param distanceThreshold The distance up to which to return an atom index.
 * @param upperTriangle Flags the inclusion of only the atoms j with j > targetAtomIndex.
 * @param sameAtomIncluded Flags the inclusion of the same atom in the set.
 * @return int The indices of the atoms within a certain distance to a given target.
 */
std::vector<int> getIndicesCloseToAtom(const PositionCollection& positions, int targetAtomIndex, double distanceThreshold,
                                       bool sameAtomIncluded = false, bool upperTriangle = true);

/**
 * @brief Builds for each atom a list of atoms within a certain distance to it.
 * @param positions A set of positions to be traversed.
 * @param distanceThreshold The distance up to which to return an atom index.
 * @param sameAtomIncluded Flags the inclusion of the same atom in the set.
 * @param upperTriangle calculate only the upper triangle of the atom pair list, i.e. j > i.
 * @return std::map<int, std::vector<int>> For each atom, the indices of the atoms within a certain distance to it
 * (upper-triangle).
 */
std::map<int, std::vector<int>> constructAtomPairList(const PositionCollection& positions, double distanceThreshold,
                                                      bool sameAtomIncluded = false, bool upperTriangle = true);

/**
 * @brief Get the index of a given atom in a given molecular structure.
 * @param structure The structure to be traversed.
 * @param atom The target atom that is tried to be found in the given structure.
 * @param distanceConsideredZero Squared distance between two atoms resulting in them being considered equal.
 * @return int The index of the given atom in the given molecular structure.
 * @throws std::runtime_error Function throws error if the atom is not found in the given structure.
 */
int getIndexOfAtomInStructure(const AtomCollection& structure, const Atom& atom, double squaredDistanceConsideredZero = 1e-4);

/**
 * @brief Get list of atoms with an absolute deviation higher than a threshold between two structures.
 * @param reference The reference positions.
 * @param positions The positions to be aligned, will be transformed in place.
 * @param elements The elements composing the structure to allow a mass-weighting of the alignment.
 * @param threshold The threshold beyond which two atoms are considered to be diverging.
 * @return A vector containing the indices of the diverging atoms.
 */
std::vector<int> getListOfDivergingAtoms(const PositionCollection& reference, PositionCollection positions,
                                         double threshold, const ElementTypeCollection& elements = {});

/**
 * @brief Get list of atoms with an absolute deviation higher than a threshold between two structures
 *        with an iterative robust alignment ignoring most flexible parts.
 * @param reference The reference positions.
 * @param positions The positions to be aligned, will be transformed in place.
 * @param threshold The threshold beyond which two atoms are considered to be diverging.
 * @param stopCriterion Threshold for the rmsd deviation of the difference of the rmsd of the alignment
 *                      between two iterations. If lower, iterative algorithm stops.
 * @param maxIterations The maximal number of iterations in the iterative alignment.
 * @param elements The elements composing the structure to allow a mass-weighting of the alignment.
 * @return A vector containing the indices of the robustly aligned diverging atoms.
 */
std::vector<int> getListOfDivergingAtomsRobust(const PositionCollection& reference, PositionCollection positions,
                                               double threshold, double stopCriterion = 1e-2, int maxIterations = 20,
                                               const ElementTypeCollection& elements = {},
                                               Core::Log log = Core::Log::silent());
} /* namespace Distances */
} // namespace Geometry
} // namespace Utils
} // namespace Scine

#endif // UTILS_GEOMETRY_DISTANCES_H_
