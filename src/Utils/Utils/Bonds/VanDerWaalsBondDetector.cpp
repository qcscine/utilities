/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Bonds/VanDerWaalsBondDetector.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include <cassert>

namespace Scine {
namespace Utils {

BondOrderCollection VanDerWaalsBondDetector::detectBonds(const AtomCollection& atoms) {
  return detectBonds(atoms.getElements(), atoms.getPositions());
}

BondOrderCollection VanDerWaalsBondDetector::detectBonds(const ElementTypeCollection& elements,
                                                         const PositionCollection& positions) {
  const int N = elements.size();
  assert(N == positions.size());
  BondOrderCollection bc(elements.size());

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < i; ++j) {
      if (vdwBondExists(elements[i], elements[j], positions.row(i), positions.row(j))) {
        bc.setOrder(i, j, 1.0);
      }
    }
  }

  return bc;
}

bool VanDerWaalsBondDetector::vdwBondExists(const ElementType& e1, const ElementType& e2, const Position& p1,
                                            const Position& p2) {
  try {
    auto vdw1 = ElementInfo::vdwRadius(e1);
    auto vdw2 = ElementInfo::vdwRadius(e2);
    auto vdwRadiiMean = (vdw1 + vdw2) / 2;
    return ((p2 - p1).squaredNorm() < vdwRadiiMean * vdwRadiiMean);
  }
  catch (...) {
    return false;
  }
}

} /* namespace Utils */
} /* namespace Scine */
