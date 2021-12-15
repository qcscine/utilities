/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Gtf.h"
#include <cmath>
#include <stdexcept>

namespace Scine {
namespace Utils {

namespace {
const double pi = 4.0 * std::atan(1.0);
}

Gtf::Gtf(int l, double expo, double coef) : exponent(expo), coefficient(coef) {
  setNormalized(l);
}

void Gtf::setNormalized(int l) {
  if (l == 0) {
    normalizedCoefficient = coefficient * std::pow(2 * exponent / pi, 0.75);
  }
  else if (l == 1) {
    normalizedCoefficient = coefficient * std::pow(2, 1.75) * std::pow(exponent, 1.25) / std::pow(pi, 0.75);
  }
  else if (l == 2) {
    normalizedCoefficient = coefficient * std::pow(2, 2.75) * std::pow(exponent, 1.75) / std::pow(pi, 0.75); // for
                                                                                                             // xy,yz,xz
                                                                                                             // orbitals
    // normalizedCoefficient2_ =
    // coefficient*std::pow(2,2.75)*std::pow(exponent,1.75)/(std::pow(pi,0.75)*std::sqrt(3.0));//for z2, r2-3z2
  }
}

} // namespace Utils
} // namespace Scine
