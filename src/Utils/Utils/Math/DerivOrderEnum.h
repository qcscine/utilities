/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_DERIVORDERENUM_H
#define UTILS_MATH_DERIVORDERENUM_H

namespace Scine {
namespace Utils {

/**
 * @class DerivativeOrder DerivOrderEnum.h
 * @brief Enum for the order up to which derivatives need to be calculated for a value.
 */
enum class DerivativeOrder { Zero, One, Two };

/**
 * @class Derivative DerivOrderEnum.h
 * @brief Enum for the derivatives to be calculated for molecular systems.
 */
enum class Derivative { None, First, SecondAtomic, SecondFull };

} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_DERIVORDERENUM_H
