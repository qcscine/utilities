/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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
template<Derivative O>
struct Derivative3DImpl { using ValueType = void; }; //Default
template<>
struct Derivative3DImpl<Derivative::First> { using ValueType = Gradient; };
template<>
struct Derivative3DImpl<Derivative::SecondAtomic> { using ValueType = Second3D; };
template<>
struct Derivative3DImpl<Derivative::SecondFull> { using ValueType = Second3D; };

// Derivative container type
template<Derivative O>
struct DerivativeContainer3DImpl { using ValueType = void; }; //Default
template<>
struct DerivativeContainer3DImpl<Derivative::First> { using ValueType = GradientCollection; };
template<>
struct DerivativeContainer3DImpl<Derivative::SecondAtomic> { using ValueType = AtomicSecondDerivativeCollection; };
template<>
struct DerivativeContainer3DImpl<Derivative::SecondFull> { using ValueType = FullSecondDerivativeCollection; };

// Underlying order of derivatives
template<Derivative O>
struct UnderlyingOrderImpl { static const DerivativeOrder o = DerivativeOrder::Zero; }; //Default
template<>
struct UnderlyingOrderImpl<Derivative::First> { static const DerivativeOrder o = DerivativeOrder::One; };
template<>
struct UnderlyingOrderImpl<Derivative::SecondAtomic> { static const DerivativeOrder o = DerivativeOrder::Two; };
template<>
struct UnderlyingOrderImpl<Derivative::SecondFull> { static const DerivativeOrder o = DerivativeOrder::Two; };
/*! \endcond */

/**
 * @brief Templated type for a derivative.
 */
template<Derivative o> using DerivativeType = typename Derivative3DImpl<o>::ValueType;
/**
 * @brief Templated type for a derivative container.
 */
template<Derivative o> using DerivativeContainerType = typename DerivativeContainer3DImpl<o>::ValueType;

/**
 * @brief Template variable for the derivative order corresponding to a derivative type.
 */
template<Derivative O>
constexpr DerivativeOrder UnderlyingOrder = UnderlyingOrderImpl<O>::o;

// clang-format on

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine

#endif // AUTOMATICDIFFERENTIATION_METHODSTYPESHELPER_H
