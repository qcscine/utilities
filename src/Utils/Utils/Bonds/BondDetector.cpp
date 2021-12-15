/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Bonds/BondDetector.h"
#include "Utils/Bonds/BondDetectorRadii.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Constants.h"
#include "Utils/Geometry/AtomCollection.h"

namespace Scine {
namespace Utils {

BondOrderCollection BondDetector::detectBonds(const AtomCollection& atoms) {
  return detectBonds(atoms.getElements(), atoms.getPositions());
}

BondOrderCollection BondDetector::detectBonds(const PeriodicSystem& periodicSystem, bool bondsAcrossBoundariesNegative) {
  return detectBonds(periodicSystem.atoms, periodicSystem.pbc, bondsAcrossBoundariesNegative);
}

BondOrderCollection BondDetector::detectBonds(const AtomCollection& atoms, const PeriodicBoundaries& pbc,
                                              bool bondsAcrossBoundariesNegative) {
  return detectBonds(atoms.getElements(), atoms.getPositions(), pbc, bondsAcrossBoundariesNegative);
}

BondOrderCollection BondDetector::detectBonds(const ElementTypeCollection& elements, const PositionCollection& positions) {
  const int N = elements.size();
  assert(N == positions.rows());
  BondOrderCollection bc(elements.size());

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < i; ++j) {
      if (bondExists(elements[i], elements[j], positions.row(i), positions.row(j))) {
        bc.setOrder(i, j, 1.0);
      }
    }
  }

  return bc;
}

BondOrderCollection BondDetector::detectBonds(const ElementTypeCollection& elements, const PositionCollection& positions,
                                              const PeriodicBoundaries& pbc, bool bondsAcrossBoundariesNegative) {
  assert(static_cast<Eigen::Index>(elements.size()) == positions.rows());
  BondOrderCollection bc(elements.size());

  const int N = elements.size();
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < i; ++j) {
      if (bondExists(elements[i], elements[j], positions.row(i), positions.row(j), pbc)) {
        double order = 1.0;
        if (bondsAcrossBoundariesNegative && pbc.minimumDistanceViaImage(positions.row(i), positions.row(j))) {
          order = -1.0;
        }
        bc.setOrder(i, j, order);
      }
    }
  }

  return bc;
}

bool BondDetector::bondExists(const ElementType& e1, const ElementType& e2, const Position& p1, const Position& p2) {
  constexpr double toleranceDistance = toBohr(Angstrom(0.4));
  auto covalentRadius1 = getCovalentRadius(e1);
  auto covalentRadius2 = getCovalentRadius(e2);

  double threshold = covalentRadius1 + covalentRadius2 + toleranceDistance;
  double thresholdSquared = threshold * threshold;

  double nucleiDistanceSquared = (p1 - p2).squaredNorm();
  return nucleiDistanceSquared < thresholdSquared;
}

bool BondDetector::bondExists(const ElementType& e1, const ElementType& e2, const Position& p1, const Position& p2,
                              const PeriodicBoundaries& pbc) {
  constexpr double toleranceDistance = toBohr(Angstrom(0.4));
  auto covalentRadius1 = getCovalentRadius(e1);
  auto covalentRadius2 = getCovalentRadius(e2);

  double threshold = covalentRadius1 + covalentRadius2 + toleranceDistance;
  double thresholdSquared = threshold * threshold;

  double nucleiDistanceSquared = Geometry::Distances::distanceSquared(p1, p2, pbc);
  return nucleiDistanceSquared < thresholdSquared;
}

double BondDetector::getCovalentRadius(ElementType e) {
  static BondDetectorRadii radii;
  return radii.getRadius(e);
}

} /* namespace Utils */
} /* namespace Scine */
