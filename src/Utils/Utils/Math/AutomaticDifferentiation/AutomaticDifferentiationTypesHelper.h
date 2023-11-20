/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef AUTOMATICDIFFERENTIATION_AUTOMATICDIFFERENTIATIONTYPESHELPER_H
#define AUTOMATICDIFFERENTIATION_AUTOMATICDIFFERENTIATIONTYPESHELPER_H

#include "../DerivOrderEnum.h"

/**
 * This header contains alias definitions defining which classes to use for the different degrees of derivatives.
 */

namespace Scine {
namespace Utils {
namespace AutomaticDifferentiation {
// Forward declaration of the derivative objects.
class First1D;
class First3D;
class Second1D;
class Second3D;

/**
 * @brief Values in 3 dimensions
 */
template<DerivativeOrder O>
struct Value3DOrder {
  using ValueType = double;
}; // Default
template<>
struct Value3DOrder<DerivativeOrder::One> {
  using ValueType = First3D;
};
template<>
struct Value3DOrder<DerivativeOrder::Two> {
  using ValueType = Second3D;
};

/**
 * @brief Values in 1 dimension
 */
template<DerivativeOrder O>
struct Value1DOrder {
  using ValueType = double;
}; // Default
template<>
struct Value1DOrder<DerivativeOrder::One> {
  using ValueType = First1D;
};
template<>
struct Value1DOrder<DerivativeOrder::Two> {
  using ValueType = Second1D;
};

/**
 * @brief Templated type for a value in 1 dimension.
 */
template<DerivativeOrder o>
using Value1DType = typename Value1DOrder<o>::ValueType;
/**
 * @brief Templated type for a value in 3 dimension.
 */
template<DerivativeOrder o>
using Value3DType = typename Value3DOrder<o>::ValueType;

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine

#endif // AUTOMATICDIFFERENTIATION_AUTOMATICDIFFERENTIATIONTYPESHELPER_H
