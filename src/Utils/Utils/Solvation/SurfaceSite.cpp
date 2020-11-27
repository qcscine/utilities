/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Solvation/SurfaceSite.h"

namespace Scine {
namespace Utils {

MolecularSurface::SurfaceSite::SurfaceSite(const Position& surfPos, const Position& orgPos) {
  position = surfPos;
  // vector from surface point normal
  normal = (surfPos - orgPos).normalized();
}

} /* namespace Utils */
} /* namespace Scine */
