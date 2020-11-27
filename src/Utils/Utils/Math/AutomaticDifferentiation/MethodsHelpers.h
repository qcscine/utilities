/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef AUTOMATICDIFFERENTIATION_METHODSHELPERS_H
#define AUTOMATICDIFFERENTIATION_METHODSHELPERS_H

#include "../AtomicSecondDerivativeCollection.h"
#include "../DerivOrderEnum.h"
#include "../FullSecondDerivativeCollection.h"
#include "AutomaticDifferentiationHelpers.h"
#include "MethodsTypesHelper.h"
#include "TypeDefinitions.h"
#include <Eigen/Core>

namespace Scine {
namespace Utils {
namespace AutomaticDifferentiation {

/**
 * @brief Get opposed derivatives. i.e. sign of first derivative changes.
 */
template<Derivative o>
DerivativeType<o> getOppositeDerivative(const DerivativeType<o>& v);
/**
 * @brief Extract the derivatives from a value type with derivatives.
 */
template<Derivative o>
DerivativeType<o> getDerivativeFromValueWithDerivatives(const Value3DType<UnderlyingOrder<o>>& v);

/**
 * @brief Add derivatives to a derivative container.
 * @param v Derivative for the pair of indexes a-b. Considering the direction a->b.
 */
template<Derivative o>
void addDerivativeToContainer(DerivativeContainerType<o>& container, int a, int b, const DerivativeType<o>& v);

/*
 * Inline implementations of all functions declared above.
 */

template<>
inline DerivativeType<Derivative::First> getOppositeDerivative<Derivative::First>(const DerivativeType<Derivative::First>& v) {
  return -v;
}
template<>
inline DerivativeType<Derivative::SecondAtomic>
getOppositeDerivative<Derivative::SecondAtomic>(const DerivativeType<Derivative::SecondAtomic>& v) {
  return v.opposite();
}
template<>
inline DerivativeType<Derivative::SecondFull>
getOppositeDerivative<Derivative::SecondFull>(const DerivativeType<Derivative::SecondFull>& v) {
  return v.opposite();
}

template<>
inline DerivativeType<Derivative::First>
getDerivativeFromValueWithDerivatives<Derivative::First>(const Value3DType<DerivativeOrder::One>& v) {
  return Gradient(v.derivatives());
}
template<>
inline DerivativeType<Derivative::SecondAtomic>
getDerivativeFromValueWithDerivatives<Derivative::SecondAtomic>(const Value3DType<DerivativeOrder::Two>& v) {
  return v;
}
template<>
inline DerivativeType<Derivative::SecondFull>
getDerivativeFromValueWithDerivatives<Derivative::SecondFull>(const Value3DType<DerivativeOrder::Two>& v) {
  return v;
}

template<>
inline void addDerivativeToContainer<Derivative::First>(DerivativeContainerType<Derivative::First>& container, int a,
                                                        int b, const DerivativeType<Derivative::First>& v) {
  container.row(b) += v;
  container.row(a) -= v;
}
template<>
inline void addDerivativeToContainer<Derivative::SecondAtomic>(DerivativeContainerType<Derivative::SecondAtomic>& container,
                                                               int a, int b,
                                                               const DerivativeType<Derivative::SecondAtomic>& v) {
  container[b] += v;
  container[a] += getOppositeDerivative<Derivative::SecondAtomic>(v);
}
template<>
inline void addDerivativeToContainer<Derivative::SecondFull>(DerivativeContainerType<Derivative::SecondFull>& container,
                                                             int a, int b, const DerivativeType<Derivative::SecondFull>& v) {
  container.addDerivative(a, b, v);
}

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine

#endif // AUTOMATICDIFFERENTIATION_METHODSHELPERS_H
