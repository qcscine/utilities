/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Second3D.h"

namespace Scine {
namespace Utils {
namespace AutomaticDifferentiation {

// Implementation of the differentiation of square root functions
Second3D sqrt(const Second3D& der) {
  double sq = std::sqrt(der.value());
  double inv2sq = 0.5 / sq;
  double inv4sq3 = 2 * inv2sq * inv2sq * inv2sq;
  return {sq,
          inv2sq * der.dx(),
          inv2sq * der.dy(),
          inv2sq * der.dz(),
          (2 * der.value() * der.XX() - der.dx() * der.dx()) * inv4sq3,
          (2 * der.value() * der.YY() - der.dy() * der.dy()) * inv4sq3,
          (2 * der.value() * der.ZZ() - der.dz() * der.dz()) * inv4sq3,
          (2 * der.value() * der.XY() - der.dx() * der.dy()) * inv4sq3,
          (2 * der.value() * der.XZ() - der.dx() * der.dz()) * inv4sq3,
          (2 * der.value() * der.YZ() - der.dy() * der.dz()) * inv4sq3};
}

// Implementation of the differentiation of exponential functions
Second3D exp(const Second3D& der) {
  double ex = std::exp(der.value());
  return {ex,
          ex * der.dx(),
          ex * der.dy(),
          ex * der.dz(),
          ex * (der.dx() * der.dx() + der.XX()),
          ex * (der.dy() * der.dy() + der.YY()),
          ex * (der.dz() * der.dz() + der.ZZ()),
          ex * (der.dx() * der.dy() + der.XY()),
          ex * (der.dx() * der.dz() + der.XZ()),
          ex * (der.dy() * der.dz() + der.YZ())};
}

// Implementation of the differentiation of the arccos function
Second3D arccos(const Second3D& der) {
  double square = der.value() * der.value();
  double root = std::sqrt(1 - square);
  double invroot = 1.0 / root;
  double invroot3 = invroot * invroot * invroot;

  return {std::acos(der.value()),
          -der.dx() * invroot,
          -der.dy() * invroot,
          -der.dz() * invroot,
          (square * der.XX() - der.XX() - der.value() * der.dx() * der.dx()) * invroot3,
          (square * der.YY() - der.YY() - der.value() * der.dy() * der.dy()) * invroot3,
          (square * der.ZZ() - der.ZZ() - der.value() * der.dz() * der.dz()) * invroot3,
          (square * der.XY() - der.XY() - der.value() * der.dx() * der.dy()) * invroot3,
          (square * der.XZ() - der.XZ() - der.value() * der.dx() * der.dz()) * invroot3,
          (square * der.YZ() - der.YZ() - der.value() * der.dy() * der.dz()) * invroot3};
}

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine
