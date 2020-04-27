/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef AUTOMATICDIFFERENTIATION_FIRSTBASE_H
#define AUTOMATICDIFFERENTIATION_FIRSTBASE_H

#include <Eigen/Core>
#include <cmath>

namespace Scine {
namespace Utils {
namespace AutomaticDifferentiation {

/**
 * @class FirstBase FirstBase.h
 * @brief Base class representing values with some derivative type and allowing for the automatic calculation of first
 *        derivatives.
 *
 *  This class makes use of CRTP so that the functions will directly have the correct return type for the
 *  derived classes.
 *
 */
template<typename DerivativeT, typename Crtp>
class FirstBase {
 public:
  /**
   * @brief Set the value of the function.
   */
  void setValue(double v);
  /**
   * @brief Get the value of the function.
   * @return double
   */
  double value() const;
  /**
   * @brief Get derivatives.
   * @return DerivativeT& Object type of derivatives in the derived class.
   */
  const DerivativeT& derivatives() const;

  /**
   * @brief Definition of operators for differentiation rules.
   */
  const Crtp& operator+=(double v);
  const Crtp& operator-=(double v);
  const Crtp& operator+=(const FirstBase& rhs);
  const Crtp& operator-=(const FirstBase& rhs);
  const Crtp& operator*=(const FirstBase& rhs);
  const Crtp& operator/=(const FirstBase& rhs);
  const Crtp& operator*=(double rhs);
  const Crtp& operator/=(double rhs);

  Crtp operator+() const;
  Crtp operator-() const;

  Crtp operator+(double f) const;
  Crtp operator-(double f) const;
  Crtp operator+(const FirstBase& rhs) const;
  Crtp operator-(const FirstBase& rhs) const;
  Crtp operator*(double f) const;
  Crtp operator/(double f) const;
  Crtp operator*(const FirstBase& rhs) const;
  Crtp operator/(const FirstBase& rhs) const;

  /**
   * @brief Returns the opposite of the derivative.
   * @return Crtp Type of the derived class
   */
  Crtp opposite() const;

 protected:
  FirstBase(double v, DerivativeT d);

  /* value v_ and derivative d_ of the function */
  double v_;
  DerivativeT d_;

 private:
  /* corresponds to "*this", but returns the derived type. */
  Crtp& thisImpl();
  /* corresponds to "*this", but returns the derived type. */
  const Crtp& thisImpl() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename DerivativeT, typename Crtp>
Crtp operator+(double v, const FirstBase<DerivativeT, Crtp>& rhs);
template<typename DerivativeT, typename Crtp>
Crtp operator-(double v, const FirstBase<DerivativeT, Crtp>& rhs);
template<typename DerivativeT, typename Crtp>
Crtp operator*(double f, const FirstBase<DerivativeT, Crtp>& rhs);
template<typename DerivativeT, typename Crtp>
Crtp operator/(double f, const FirstBase<DerivativeT, Crtp>& rhs);

template<typename DerivativeT, typename Crtp>
Crtp square(const FirstBase<DerivativeT, Crtp>& value);
template<typename DerivativeT, typename Crtp>
Crtp sqrt(const FirstBase<DerivativeT, Crtp>& value);
template<typename DerivativeT, typename Crtp>
Crtp exp(const FirstBase<DerivativeT, Crtp>& value);

template<typename DerivativeT, typename Crtp>
FirstBase<DerivativeT, Crtp>::FirstBase(double v, DerivativeT d) : v_(v), d_(std::move(d)) {
}

template<typename DerivativeT, typename Crtp>
Crtp& FirstBase<DerivativeT, Crtp>::thisImpl() {
  return static_cast<Crtp&>(*this);
}

template<typename DerivativeT, typename Crtp>
const Crtp& FirstBase<DerivativeT, Crtp>::thisImpl() const {
  return static_cast<const Crtp&>(*this);
}

/*
 *
 * Implementation of the differentiation rules as inline functions
 *
 */

template<typename DerivativeT, typename Crtp>
void FirstBase<DerivativeT, Crtp>::setValue(double v) {
  v_ = v;
}

template<typename DerivativeT, typename Crtp>
const DerivativeT& FirstBase<DerivativeT, Crtp>::derivatives() const {
  return d_;
}

template<typename DerivativeT, typename Crtp>
double FirstBase<DerivativeT, Crtp>::value() const {
  return v_;
}

template<typename DerivativeT, typename Crtp>
Crtp FirstBase<DerivativeT, Crtp>::operator+() const {
  return {v_, d_};
}

template<typename DerivativeT, typename Crtp>
Crtp FirstBase<DerivativeT, Crtp>::operator-() const {
  return {-v_, -d_};
}

template<typename DerivativeT, typename Crtp>
const Crtp& FirstBase<DerivativeT, Crtp>::operator+=(double v) {
  v_ += v;
  return thisImpl();
}

template<typename DerivativeT, typename Crtp>
const Crtp& FirstBase<DerivativeT, Crtp>::operator-=(double v) {
  v_ -= v;
  return thisImpl();
}

template<typename DerivativeT, typename Crtp>
const Crtp& FirstBase<DerivativeT, Crtp>::operator+=(const FirstBase<DerivativeT, Crtp>& rhs) {
  v_ += rhs.v_;
  d_ += rhs.d_;
  return thisImpl();
}

template<typename DerivativeT, typename Crtp>
const Crtp& FirstBase<DerivativeT, Crtp>::operator-=(const FirstBase<DerivativeT, Crtp>& rhs) {
  v_ -= rhs.v_;
  d_ -= rhs.d_;
  return thisImpl();
}

template<typename DerivativeT, typename Crtp>
const Crtp& FirstBase<DerivativeT, Crtp>::operator*=(const FirstBase<DerivativeT, Crtp>& rhs) {
  d_ = v_ * rhs.d_ + rhs.v_ * d_;
  v_ *= rhs.v_;
  return thisImpl();
}

template<typename DerivativeT, typename Crtp>
const Crtp& FirstBase<DerivativeT, Crtp>::operator/=(const FirstBase<DerivativeT, Crtp>& rhs) {
  d_ = d_ / rhs.v_ - v_ / (rhs.v_ * rhs.v_) * rhs.d_;
  v_ /= rhs.v_;
  return thisImpl();
}

template<typename DerivativeT, typename Crtp>
const Crtp& FirstBase<DerivativeT, Crtp>::operator*=(double rhs) {
  v_ *= rhs;
  d_ *= rhs;
  return thisImpl();
}

template<typename DerivativeT, typename Crtp>
const Crtp& FirstBase<DerivativeT, Crtp>::operator/=(double rhs) {
  v_ /= rhs;
  d_ /= rhs;
  return thisImpl();
}

template<typename DerivativeT, typename Crtp>
Crtp FirstBase<DerivativeT, Crtp>::operator+(double f) const {
  auto a = thisImpl();
  a += f;
  return a;
}

template<typename DerivativeT, typename Crtp>
Crtp FirstBase<DerivativeT, Crtp>::operator-(double f) const {
  auto a = thisImpl();
  a -= f;
  return a;
}

template<typename DerivativeT, typename Crtp>
Crtp FirstBase<DerivativeT, Crtp>::operator+(const FirstBase<DerivativeT, Crtp>& rhs) const {
  auto f = thisImpl();
  f += rhs;
  return f;
}

template<typename DerivativeT, typename Crtp>
Crtp FirstBase<DerivativeT, Crtp>::operator-(const FirstBase<DerivativeT, Crtp>& rhs) const {
  auto f = thisImpl();
  f -= rhs;
  return f;
}

template<typename DerivativeT, typename Crtp>
Crtp FirstBase<DerivativeT, Crtp>::operator*(double f) const {
  auto a = thisImpl();
  a *= f;
  return a;
}

template<typename DerivativeT, typename Crtp>
Crtp FirstBase<DerivativeT, Crtp>::operator*(const FirstBase<DerivativeT, Crtp>& rhs) const {
  auto a = thisImpl();
  a *= rhs;
  return a;
}

template<typename DerivativeT, typename Crtp>
Crtp FirstBase<DerivativeT, Crtp>::operator/(double f) const {
  auto a = thisImpl();
  a /= f;
  return a;
}

template<typename DerivativeT, typename Crtp>
Crtp FirstBase<DerivativeT, Crtp>::operator/(const FirstBase<DerivativeT, Crtp>& rhs) const {
  auto a = thisImpl();
  a /= rhs;
  return a;
}

template<typename DerivativeT, typename Crtp>
Crtp FirstBase<DerivativeT, Crtp>::opposite() const {
  return {v_, -d_};
}

template<typename DerivativeT, typename Crtp>
Crtp operator+(double v, const FirstBase<DerivativeT, Crtp>& rhs) {
  return rhs + v;
}

template<typename DerivativeT, typename Crtp>
Crtp operator-(double v, const FirstBase<DerivativeT, Crtp>& rhs) {
  auto val = -rhs;
  val += v;
  return val;
}

template<typename DerivativeT, typename Crtp>
Crtp operator*(double f, const FirstBase<DerivativeT, Crtp>& rhs) {
  return rhs * f;
}

template<typename DerivativeT, typename Crtp>
Crtp operator/(double f, const FirstBase<DerivativeT, Crtp>& rhs) {
  return {f / rhs.value(), -f / (rhs.value() * rhs.value()) * rhs.derivatives()};
}

// deriving the square root function
template<typename DerivativeT, typename Crtp>
Crtp sqrt(const FirstBase<DerivativeT, Crtp>& value) {
  double sq = std::sqrt(value.value());
  double invsq = 0.5 / sq;
  return {sq, invsq * value.derivatives()};
}

// deriving the square function
template<typename DerivativeT, typename Crtp>
Crtp square(const FirstBase<DerivativeT, Crtp>& value) {
  double square = value.value() * value.value();
  return {square, 2 * value.value() * value.derivatives()};
}

// deriving the exponential function
template<typename DerivativeT, typename Crtp>
Crtp exp(const FirstBase<DerivativeT, Crtp>& value) {
  double e = std::exp(value.value());
  return {e, e * value.derivatives()};
}

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine

#endif // AUTOMATICDIFFERENTIATION_FIRSTBASE_H
