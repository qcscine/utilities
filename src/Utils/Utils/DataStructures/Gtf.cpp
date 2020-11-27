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
const double pi = 4.0 * atan(1.0);
}

Gtf::Gtf(int l, double expo, double coef) : exponent(expo), coefficient(coef) {
  setNormalized(l);
}

void Gtf::setNormalized(int l) {
  if (l == 0) {
    normalizedCoefficient = coefficient * pow(2 * exponent / pi, 0.75);
  }
  else if (l == 1) {
    normalizedCoefficient = coefficient * pow(2, 1.75) * pow(exponent, 1.25) / pow(pi, 0.75);
  }
  else if (l == 2) {
    normalizedCoefficient = coefficient * pow(2, 2.75) * pow(exponent, 1.75) / pow(pi, 0.75); // for xy,yz,xz orbitals
    // normalizedCoefficient2_ = coefficient*pow(2,2.75)*pow(exponent,1.75)/(pow(pi,0.75)*sqrt(3.0));//for z2, r2-3z2
  }
}

} // namespace Utils
} // namespace Scine
