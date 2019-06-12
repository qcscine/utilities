/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_DERIVORDERENUM_H
#define UTILS_MATH_DERIVORDERENUM_H

namespace Scine {
namespace Utils {

/**
 * @class derivOrder DerivOrderEnum.h
 * @brief Enum for the order up to which derivatives need to be calculated for a value.
 */
enum class derivOrder { zero, one, two };

/**
 * @class derivativeType DerivOrderEnum.h
 * @brief Enum for the derivatives to be calculated for molecular systems.
 */
enum class derivativeType { none, first, second_atomic, second_full };

} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_DERIVORDERENUM_H
