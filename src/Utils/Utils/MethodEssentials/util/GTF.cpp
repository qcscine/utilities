/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GTF.h"
#include <cmath>

namespace Scine {
namespace Utils {

namespace {
const double pi = 4.0 * atan(1.0);
}

void GTF::set(int l, double a, double c) {
  exponent_ = a;
  coefficient_ = c;

  normalizeCoefficient(l, c);
}

void GTF::normalizeCoefficient(int l, double c) {
  if (l == 0) {
    normalizedCoefficient_ = c * pow(2 * exponent_ / pi, 0.75);
  }
  else if (l == 1) {
    normalizedCoefficient_ = c * pow(2, 1.75) * pow(exponent_, 1.25) / pow(pi, 0.75);
  }
  else if (l == 2) {
    normalizedCoefficient_ = c * pow(2, 2.75) * pow(exponent_, 1.75) / pow(pi, 0.75); // for xy,yz,xz orbitals
    // normalizedCoefficient2_ = c*pow(2,2.75)*pow(exponent_,1.75)/(pow(pi,0.75)*sqrt(3.0));//for z2, r2-3z2
  }
}

} // namespace Utils
} // namespace Scine
