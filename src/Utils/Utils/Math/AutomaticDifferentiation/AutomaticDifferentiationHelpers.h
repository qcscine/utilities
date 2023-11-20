/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef AUTOMATICDIFFERENTIATION_AUTOMATICDIFFERENTIATIONHELPERS_H
#define AUTOMATICDIFFERENTIATION_AUTOMATICDIFFERENTIATIONHELPERS_H

#include "../DerivOrderEnum.h"
#include "AutomaticDifferentiationTypesHelper.h"
#include "TypeDefinitions.h"
#include <Eigen/Core>

/**
 *
 * This header file contains functions that allow for common notation for common
 * things that can be done at a different degree of derivatives.
 *
 */

namespace Scine {
namespace Utils {
namespace AutomaticDifferentiation {

/**
 * @brief Get the enum corresponding to the derivative order.
 */
inline DerivativeOrder getDerivativeOrderEnum(unsigned order) {
  return (order == 0) ? DerivativeOrder::Zero : order == 1 ? DerivativeOrder::One : DerivativeOrder::Two;
}

/**
 * @brief Get the integer corresponding to the derivative order.
 */
inline int getDerivativeOrder(DerivativeOrder order) {
  return (order == DerivativeOrder::Zero) ? 0 : order == DerivativeOrder::One ? 1 : 2;
}

/**
 * @brief Create a value with derivatives in 1 dimension and neglect the unnecessary derivatives.
 */
template<DerivativeOrder o>
Value1DType<o> getFromFull(double v, double firstDer, double secondDer);
/**
 * @brief Transform a double to a ValueWithDerivative in one dimension, with derivatives equal to zero.
 */
template<DerivativeOrder o>
Value1DType<o> constant1D(double c);
/**
 * @brief Transform v to a ValueWithDerivative in one dimension, with first derivative 1 and second derivative 0.
 */
template<DerivativeOrder o>
Value1DType<o> variableWithUnitDerivative(double v);
/**
 * @brief Transform a double to a ValueWithDerivative in three dimensions, with derivatives equal to zero.
 */
template<DerivativeOrder o>
Value3DType<o> constant3D(double c);
/**
 * @brief Get X as a value with derivatives in three dimensions. The only non-zero derivative is the first derivative in
 *        x-direction.
 */
template<DerivativeOrder o>
Value3DType<o> toX(double x);
/**
 * @brief Get Y as a value with derivatives in three dimensions. The only non-zero derivative is the first derivative in
 *        y-direction.
 */
template<DerivativeOrder o>
Value3DType<o> toY(double y);
/** @brief Get Z as a value with derivatives in three dimensions. The only non-zero derivative is the first derivative
 * in z-direction.
 */
template<DerivativeOrder o>
Value3DType<o> toZ(double z);
/**
 * @brief Get R-squared as a value with derivatives in three dimensions.
 */
template<DerivativeOrder o>
Value3DType<o> toRSquared(double x, double y, double z);
/**
 * @brief Get a value with derivatives in 3 dimensions from the value with derivatives in one dimension, given a vector
 * R.
 */
template<DerivativeOrder o>
Value3DType<o> get3Dfrom1D(Value1DType<o> v, const Eigen::Vector3d& R);
/**
 * @brief Get the value with inverted derivatives (useful for pairs of derivatives for two atoms). The sign of the first
 *        derivatives changes.
 */
template<DerivativeOrder o>
Value3DType<o> getValueWithOppositeDerivative(const Value3DType<o>& v);

/**
 * @brief Extract the value with derivatives in 1 dimension as a double.
 */
template<DerivativeOrder o>
double getValue1DAsDouble(const Value1DType<o>& v);
/**
 * @brief Extract the value with derivatives in 3 dimension as a double.
 */
template<DerivativeOrder o>
double getValue3DAsDouble(const Value3DType<o>& v);

/*
 * Inline implementations of all functions declared above.
 */

template<>
inline double constant1D<DerivativeOrder::Zero>(double c) {
  return c;
}
template<>
inline double variableWithUnitDerivative<DerivativeOrder::Zero>(double v) {
  return v;
}
template<>
inline double getFromFull<DerivativeOrder::Zero>(double v, double /*firstDer*/, double /*secondDer*/) {
  return v;
}
template<>
inline double constant3D<DerivativeOrder::Zero>(double c) {
  return c;
}
template<>
inline double toX<DerivativeOrder::Zero>(double x) {
  return x;
}
template<>
inline double toY<DerivativeOrder::Zero>(double y) {
  return y;
}
template<>
inline double toZ<DerivativeOrder::Zero>(double z) {
  return z;
}
template<>
inline double toRSquared<DerivativeOrder::Zero>(double x, double y, double z) {
  return x * x + y * y + z * z;
}
template<>
inline double get3Dfrom1D<DerivativeOrder::Zero>(double v, const Eigen::Vector3d& /*R*/) {
  return v;
}
template<>
inline double getValueWithOppositeDerivative<DerivativeOrder::Zero>(const double& v) {
  return v;
}
template<>
inline double getValue1DAsDouble<DerivativeOrder::Zero>(const double& v) {
  return v;
}
template<>
inline double getValue3DAsDouble<DerivativeOrder::Zero>(const double& v) {
  return v;
}

template<>
inline First1D constant1D<DerivativeOrder::One>(double c) {
  return {c, 0};
}
template<>
inline First1D variableWithUnitDerivative<DerivativeOrder::One>(double v) {
  return {v, 1};
}
template<>
inline First1D getFromFull<DerivativeOrder::One>(double v, double firstDer, double /*secondDer*/) {
  return {v, firstDer};
}
template<>
inline First3D constant3D<DerivativeOrder::One>(double c) {
  return {c, 0, 0, 0};
}
template<>
inline First3D toX<DerivativeOrder::One>(double x) {
  return {x, 1, 0, 0};
}
template<>
inline First3D toY<DerivativeOrder::One>(double y) {
  return {y, 0, 1, 0};
}
template<>
inline First3D toZ<DerivativeOrder::One>(double z) {
  return {z, 0, 0, 1};
}
template<>
inline First3D toRSquared<DerivativeOrder::One>(double x, double y, double z) {
  return {x * x + y * y + z * z, 2 * x, 2 * y, 2 * z};
}
template<>
inline First3D get3Dfrom1D<DerivativeOrder::One>(First1D v, const Eigen::Vector3d& R) {
  Eigen::Vector3d normalizedR = R.normalized();
  return {v.value(), v.derivative() * normalizedR.x(), v.derivative() * normalizedR.y(), v.derivative() * normalizedR.z()};
}
template<>
inline First3D getValueWithOppositeDerivative<DerivativeOrder::One>(const First3D& v) {
  return v.opposite();
}
template<>
inline double getValue1DAsDouble<DerivativeOrder::One>(const First1D& v) {
  return v.value();
}
template<>
inline double getValue3DAsDouble<DerivativeOrder::One>(const First3D& v) {
  return v.value();
}

template<>
inline Second1D constant1D<DerivativeOrder::Two>(double c) {
  return {c, 0, 0};
}
template<>
inline Second1D variableWithUnitDerivative<DerivativeOrder::Two>(double v) {
  return {v, 1, 0};
}
template<>
inline Second1D getFromFull<DerivativeOrder::Two>(double v, double firstDer, double secondDer) {
  return {v, firstDer, secondDer};
}
template<>
inline Second3D constant3D<DerivativeOrder::Two>(double c) {
  return {c, 0, 0, 0};
}
template<>
inline Second3D toX<DerivativeOrder::Two>(double x) {
  return {x, 1, 0, 0};
}
template<>
inline Second3D toY<DerivativeOrder::Two>(double y) {
  return {y, 0, 1, 0};
}
template<>
inline Second3D toZ<DerivativeOrder::Two>(double z) {
  return {z, 0, 0, 1};
}
template<>
inline Second3D toRSquared<DerivativeOrder::Two>(double x, double y, double z) {
  return {x * x + y * y + z * z, 2 * x, 2 * y, 2 * z, 2, 2, 2, 0, 0, 0};
}
template<>
inline Second3D get3Dfrom1D<DerivativeOrder::Two>(Second1D v, const Eigen::Vector3d& R) {
  double norm = R.norm();
  Eigen::Vector3d normalizedR = R / norm;
  double firstDerivDividedByR = v.first() / norm;
  return {v.value(),
          v.first() * normalizedR.x(),
          v.first() * normalizedR.y(),
          v.first() * normalizedR.z(),
          v.second() * normalizedR.x() * normalizedR.x() + firstDerivDividedByR * (1 - normalizedR.x() * normalizedR.x()),
          v.second() * normalizedR.y() * normalizedR.y() + firstDerivDividedByR * (1 - normalizedR.y() * normalizedR.y()),
          v.second() * normalizedR.z() * normalizedR.z() + firstDerivDividedByR * (1 - normalizedR.z() * normalizedR.z()),
          v.second() * normalizedR.x() * normalizedR.y() - firstDerivDividedByR * normalizedR.x() * normalizedR.y(),
          v.second() * normalizedR.x() * normalizedR.z() - firstDerivDividedByR * normalizedR.x() * normalizedR.z(),
          v.second() * normalizedR.y() * normalizedR.z() - firstDerivDividedByR * normalizedR.y() * normalizedR.z()};
}
template<>
inline Second3D getValueWithOppositeDerivative<DerivativeOrder::Two>(const Second3D& v) {
  return v.opposite();
}
template<>
inline double getValue1DAsDouble<DerivativeOrder::Two>(const Second1D& v) {
  return v.value();
}
template<>
inline double getValue3DAsDouble<DerivativeOrder::Two>(const Second3D& v) {
  return v.value();
}

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine

#endif // AUTOMATICDIFFERENTIATION_AUTOMATICDIFFERENTIATIONHELPERS_H
