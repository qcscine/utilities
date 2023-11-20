/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_SURFACESITE_H
#define UTILSOS_SURFACESITE_H

#include "Utils/Typenames.h"

namespace Scine {
namespace Utils {

namespace MolecularSurface {
/**
 * @class SurfaceSite SurfaceSite.h
 * @brief A surface site having position and its normal vector
 *
 * Small class for surface site on a sphere containing the position and the normal vector,
 * calculated using the origin of the sphere.
 *
 */
class SurfaceSite {
 public:
  /**
   * @brief Construct a new surface site object
   *
   * The surface site will be placed at (0,0,1) and the origin of the sphere at (0,0,0) if no input is given.
   *
   * @param surfPos Position of the surface site.
   * @param orgPos Position of the origin of the sphere.
   */
  explicit SurfaceSite(const Position& surfPos = Position(0, 0, 1), const Position& orgPos = Position(0, 0, 0));

  Position position;
  Position normal;
};

} /* namespace MolecularSurface */
} /* namespace Utils */
} /* namespace Scine */

#endif // UTILSOS_SURFACESITE_H
