/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "SolidStateBondDetector.h"
#include "BondDetector.h"
#include "BondOrderCollection.h"
#include "Utils/Constants.h"
#include "Utils/Geometry/AtomCollection.h"

namespace Scine {
namespace Utils {

BondOrderCollection SolidStateBondDetector::detectBonds(const AtomCollection& atoms,
                                                        const std::unordered_set<unsigned>& solidStateIndices) {
  return detectBonds(atoms.getElements(), atoms.getPositions(), solidStateIndices);
}

BondOrderCollection SolidStateBondDetector::detectBonds(const PeriodicSystem& periodicSystem, bool bondsAcrossBoundariesNegative) {
  return detectBonds(periodicSystem.atoms, periodicSystem.pbc, periodicSystem.solidStateAtomIndices,
                     bondsAcrossBoundariesNegative);
}

BondOrderCollection SolidStateBondDetector::detectBonds(const AtomCollection& atoms, const PeriodicBoundaries& pbc,
                                                        const std::unordered_set<unsigned>& solidStateIndices,
                                                        bool bondsAcrossBoundariesNegative) {
  return detectBonds(atoms.getElements(), atoms.getPositions(), pbc, solidStateIndices, bondsAcrossBoundariesNegative);
}

BondOrderCollection SolidStateBondDetector::detectBonds(const ElementTypeCollection& elements,
                                                        const PositionCollection& positions,
                                                        const std::unordered_set<unsigned>& solidStateIndices) {
  assert(static_cast<Eigen::Index>(elements.size()) == positions.rows());
  const long N = positions.rows();

  BondOrderCollection nnBonds = Geometry::Distances::nearestNeighborsBondOrders(positions);
  BondOrderCollection radiiBonds = BondDetector::detectBonds(elements, positions);
  BondOrderCollection final = BondOrderCollection(static_cast<int>(N));

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < i; ++j) {
      // 2 non solid state atoms -> simply decide via radii
      if (solidStateIndices.find(i) == solidStateIndices.end() && solidStateIndices.find(j) == solidStateIndices.end()) {
        final.setOrder(i, j, radiiBonds.getOrder(i, j));
      }
      // 2 solid state atoms -> simply decide via NN
      else if (solidStateIndices.find(i) != solidStateIndices.end() && solidStateIndices.find(j) != solidStateIndices.end()) {
        final.setOrder(i, j, nnBonds.getOrder(i, j));
      }
      // mixture, use radii but reevaluate NN of solid state atom
      else {
        final.setOrder(i, j, radiiBonds.getOrder(i, j));
        // if bond via NN, construct a position without the non-solid counterpart to reevaluate NN of solid state atom
        if (nnBonds.getOrder(i, j) > 0.0) {
          // either i or j is solid state
          bool jIsSolid = solidStateIndices.find(j) != solidStateIndices.end();
          int solid = (jIsSolid) ? j : i;
          int nonSolid = (jIsSolid) ? i : j;
          PositionCollection fakePositions(N - 1, 3);
          if (N == 2) {
            fakePositions << positions.row(solid);
          }
          else {
            PositionCollection topRows = positions.topRows(nonSolid);
            PositionCollection bottomRows = positions.bottomRows(N - nonSolid - 1);
            fakePositions << topRows, bottomRows;
          }
          Position solidStatePos = positions.row(solid);
          std::vector<int> neighbors = Geometry::Distances::nearestNeighborsInPositions(fakePositions, solidStatePos);
          // set order to one for neighbors if they are also solid state atoms
          for (auto neighborIndex : neighbors) {
            int realIndex = (neighborIndex < nonSolid) ? neighborIndex : neighborIndex + 1;
            if (solidStateIndices.find(static_cast<unsigned>(realIndex)) != solidStateIndices.end()) {
              assert(solid != realIndex);
              final.setOrder(solid, realIndex, 1.0);
            }
          }
        }
      }
    }
  }

  return final;
}

BondOrderCollection SolidStateBondDetector::detectBonds(const ElementTypeCollection& elements,
                                                        const PositionCollection& positions, const PeriodicBoundaries& pbc,
                                                        const std::unordered_set<unsigned>& solidStateIndices,
                                                        bool bondsAcrossBoundariesNegative) {
  assert(static_cast<Eigen::Index>(elements.size()) == positions.rows());
  const long N = positions.rows();

  BondOrderCollection nnBonds = Geometry::Distances::nearestNeighborsBondOrders(positions, pbc);
  BondOrderCollection radiiBonds = BondDetector::detectBonds(elements, positions, pbc, bondsAcrossBoundariesNegative);
  BondOrderCollection final = BondOrderCollection(static_cast<int>(N));

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < i; ++j) {
      // 2 non solid state atoms -> simply decide via radii
      if (solidStateIndices.find(i) == solidStateIndices.end() && solidStateIndices.find(j) == solidStateIndices.end()) {
        final.setOrder(i, j, radiiBonds.getOrder(i, j));
      }
      // 2 solid state atoms -> simply decide via NN, but we have to adjust possible negative bond order
      else if (solidStateIndices.find(i) != solidStateIndices.end() && solidStateIndices.find(j) != solidStateIndices.end()) {
        if (nnBonds.getOrder(i, j) > 0.0) {
          double order =
              (bondsAcrossBoundariesNegative && pbc.minimumDistanceViaImage(positions.row(i), positions.row(j))) ? -1.0 : 1.0;
          final.setOrder(i, j, order);
        }
        else {
          final.setOrder(i, j, 0.0);
        }
      }
      // mixture, use radii but reevaluate NN of solid state atom
      else {
        final.setOrder(i, j, radiiBonds.getOrder(i, j));
        // if bond via NN, construct a position without the non-solid counterpart to reevaluate NN of solid state atom
        if (nnBonds.getOrder(i, j) > 0.0) {
          // either i or j is solid state
          bool jIsSolid = solidStateIndices.find(j) != solidStateIndices.end();
          int solid = (jIsSolid) ? j : i;
          int nonSolid = (jIsSolid) ? i : j;
          PositionCollection fakePositions(N - 1, 3);
          if (N == 2) {
            fakePositions << positions.row(solid);
          }
          else {
            PositionCollection topRows = positions.topRows(nonSolid);
            PositionCollection bottomRows = positions.bottomRows(N - nonSolid - 1);
            fakePositions << topRows, bottomRows;
          }
          Position solidStatePos = positions.row(solid);
          std::vector<int> neighbors = Geometry::Distances::nearestNeighborsInPositions(fakePositions, solidStatePos, pbc);
          // set order to one for neighbors if they are also solid state atoms
          for (auto neighborIndex : neighbors) {
            int realIndex = (neighborIndex < nonSolid) ? neighborIndex : neighborIndex + 1;
            if (solidStateIndices.find(static_cast<unsigned>(realIndex)) != solidStateIndices.end()) {
              assert(solid != realIndex);
              double order = (bondsAcrossBoundariesNegative &&
                              pbc.minimumDistanceViaImage(positions.row(solid), positions.row(realIndex)))
                                 ? -1.0
                                 : 1.0;
              final.setOrder(solid, realIndex, order);
            }
          }
        }
      }
    }
  }

  return final;
}

} /* namespace Utils */
} /* namespace Scine */
