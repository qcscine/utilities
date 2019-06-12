/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef AUTOMATICDIFFERENTIATION_METHODSTYPESHELPER_H
#define AUTOMATICDIFFERENTIATION_METHODSTYPESHELPER_H

#include "../DerivOrderEnum.h"
#include "AutomaticDifferentiationTypesHelper.h"
#include <Utils/Typenames.h>

namespace Scine {
namespace Utils {

class AtomicSecondDerivativeCollection;
class FullSecondDerivativeCollection;

namespace AutomaticDifferentiation {
class First1D;
class First3D;
class Second1D;
class Second3D;

// clang-format off

/*! \cond */
// Derivative type
template<derivativeType O>
struct Derivative3DImpl { using ValueType = void; }; //Default
template<>
struct Derivative3DImpl<derivativeType::first> { using ValueType = Gradient; };
template<>
struct Derivative3DImpl<derivativeType::second_atomic> { using ValueType = Second3D; };
template<>
struct Derivative3DImpl<derivativeType::second_full> { using ValueType = Second3D; };

// Derivative container type
template<derivativeType O>
struct DerivativeContainer3DImpl { using ValueType = void; }; //Default
template<>
struct DerivativeContainer3DImpl<derivativeType::first> { using ValueType = GradientCollection; };
template<>
struct DerivativeContainer3DImpl<derivativeType::second_atomic> { using ValueType = AtomicSecondDerivativeCollection; };
template<>
struct DerivativeContainer3DImpl<derivativeType::second_full> { using ValueType = FullSecondDerivativeCollection; };

// Underlying order of derivatives
template<derivativeType O>
struct UnderlyingOrderImpl { static const derivOrder o = derivOrder::zero; }; //Default
template<>
struct UnderlyingOrderImpl<derivativeType::first> { static const derivOrder o = derivOrder::one; };
template<>
struct UnderlyingOrderImpl<derivativeType::second_atomic> { static const derivOrder o = derivOrder::two; };
template<>
struct UnderlyingOrderImpl<derivativeType::second_full> { static const derivOrder o = derivOrder::two; };
/*! \endcond */

/**
 * @brief Templated type for a derivative.
 */
template<derivativeType o> using DerivativeType = typename Derivative3DImpl<o>::ValueType;
/**
 * @brief Templated type for a derivative container.
 */
template<derivativeType o> using DerivativeContainerType = typename DerivativeContainer3DImpl<o>::ValueType;

/**
 * @brief Template variable for the derivative order corresponding to a derivative type.
 */
template<derivativeType O>
constexpr derivOrder UnderlyingOrder = UnderlyingOrderImpl<O>::o;

// clang-format on

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine

#endif // AUTOMATICDIFFERENTIATION_METHODSTYPESHELPER_H
