/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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
template<derivativeType o>
DerivativeType<o> getOppositeDerivative(const DerivativeType<o>& v);
/**
 * @brief Extract the derivatives from a value type with derivatives.
 */
template<derivativeType o>
DerivativeType<o> getDerivativeFromValueWithDerivatives(const Value3DType<UnderlyingOrder<o>>& v);

/**
 * @brief Add derivatives to a derivative container.
 * @param v Derivative for the pair of indexes a-b. Considering the direction a->b.
 */
template<derivativeType o>
void addDerivativeToContainer(DerivativeContainerType<o>& container, int a, int b, const DerivativeType<o>& v);

/*
 * Inline implementations of all functions declared above.
 */

template<>
inline DerivativeType<derivativeType::first>
getOppositeDerivative<derivativeType::first>(const DerivativeType<derivativeType::first>& v) {
  return -v;
}
template<>
inline DerivativeType<derivativeType::second_atomic>
getOppositeDerivative<derivativeType::second_atomic>(const DerivativeType<derivativeType::second_atomic>& v) {
  return v.opposite();
}
template<>
inline DerivativeType<derivativeType::second_full>
getOppositeDerivative<derivativeType::second_full>(const DerivativeType<derivativeType::second_full>& v) {
  return v.opposite();
}

template<>
inline DerivativeType<derivativeType::first>
getDerivativeFromValueWithDerivatives<derivativeType::first>(const Value3DType<derivOrder::one>& v) {
  return Gradient(v.derivatives());
}
template<>
inline DerivativeType<derivativeType::second_atomic>
getDerivativeFromValueWithDerivatives<derivativeType::second_atomic>(const Value3DType<derivOrder::two>& v) {
  return v;
}
template<>
inline DerivativeType<derivativeType::second_full>
getDerivativeFromValueWithDerivatives<derivativeType::second_full>(const Value3DType<derivOrder::two>& v) {
  return v;
}

template<>
inline void addDerivativeToContainer<derivativeType::first>(DerivativeContainerType<derivativeType::first>& container,
                                                            int a, int b, const DerivativeType<derivativeType::first>& v) {
  container.row(b) += v;
  container.row(a) -= v;
}
template<>
inline void
addDerivativeToContainer<derivativeType::second_atomic>(DerivativeContainerType<derivativeType::second_atomic>& container,
                                                        int a, int b, const DerivativeType<derivativeType::second_atomic>& v) {
  container[b] += v;
  container[a] += getOppositeDerivative<derivativeType::second_atomic>(v);
}
template<>
inline void
addDerivativeToContainer<derivativeType::second_full>(DerivativeContainerType<derivativeType::second_full>& container,
                                                      int a, int b, const DerivativeType<derivativeType::second_full>& v) {
  container.addDerivative(a, b, v);
}

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine

#endif // AUTOMATICDIFFERENTIATION_METHODSHELPERS_H
