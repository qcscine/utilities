/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Second1D.h"
#include <cmath>

namespace Scine {
namespace Utils {
namespace AutomaticDifferentiation {

// Implementation of the differentiation rule for the square root function
Second1D sqrt(const Second1D& der) {
  double sq = std::sqrt(der.value());
  return {sq, 0.5 / sq * der.first(), (2 * der.value() * der.second() - der.first() * der.first()) / (4 * sq * sq * sq)};
}

// Implementation of the differentiation rule for the exponential function
Second1D exp(const Second1D& der) {
  double ex = std::exp(der.value());
  return {ex, ex * der.first(), ex * (der.second() + der.first() * der.first())};
}

// Implementation of the differentiation rule for the cosine function
Second1D cos(const Second1D& der) {
  return {std::cos(der.value()), -std::sin(der.value()) * der.first(),
          -std::sin(der.value()) * der.second() - std::cos(der.value()) * der.first() * der.first()};
}

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine
