/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef AUTOMATICDIFFERENTIATION_SECOND1D_H
#define AUTOMATICDIFFERENTIATION_SECOND1D_H

namespace Scine {
namespace Utils {
namespace AutomaticDifferentiation {

/**
 * @class Second1D Second1D.h
 * @brief Class representing values in one dimensions and allowing for the automatic calculation of first and second
 *        derivatives.
 */
class Second1D {
 public:
  /**
   * @brief Default constructor.
   */
  Second1D() = default;
  /**
   * @brief Constructor that takes a value, first derivative and second derivative of a function as an argument.
   * @param v Value of the function.
   * @param d First derivative of the function.
   * @param h Second derivative of the function.
   */
  Second1D(double v, double d, double h);

  /**
   * @brief Getter for the value.
   * @return double
   */
  double value() const;
  /**
   * @brief Getter for the first derivative.
   * @return double
   */
  double first() const;
  /**
   * @brief Getter for the second derivative.
   * @return double
   */
  double second() const;
  /**
   * @brief Setter for the value.
   */
  void setValue(double v);
  /**
   * @brief Setter for the first derivative.
   */
  void setFirst(double d);
  /**
   * @brief Setter for the second derivative.
   */
  void setSecond(double h);
  /**
   * @brief Zero initializer.
   */
  void setZero();

  /**
   *  Definition of all differentiation rules.
   */
  const Second1D& operator+=(const Second1D& rhs);
  const Second1D& operator-=(const Second1D& rhs);
  const Second1D& operator*=(const Second1D& rhs);
  const Second1D& operator/=(const Second1D& rhs);
  const Second1D& operator*=(double rhs);
  const Second1D& operator/=(double rhs);

  Second1D operator+();
  Second1D operator-();
  Second1D operator+(const Second1D& der) const;
  Second1D operator-(const Second1D& rhs) const;
  Second1D operator*(const Second1D& rhs) const;
  Second1D operator/(const Second1D& rhs) const;
  Second1D operator*(double f) const;
  Second1D operator/(double f) const;

 private:
  // value of the function
  double v_{0};
  // first derivative of the function
  double d_{0};
  // second derivative of the function
  double h_{0};
};

Second1D operator*(double f, const Second1D& der);
Second1D operator+(double f, const Second1D& h);
Second1D operator+(const Second1D& h, double f);
Second1D operator-(double f, const Second1D& h);
Second1D operator-(const Second1D& h, double f);
Second1D operator/(double v, const Second1D& der);

Second1D sqrt(const Second1D& der);
Second1D exp(const Second1D& der);
Second1D cos(const Second1D& der);

// Constructor that takes the value v, first derivative d and second derivative h as an argument.
inline Second1D::Second1D(double v, double d, double h) : v_(v), d_(d), h_(h) {
}

/*
 *
 * Inline implementation of getters and setters.
 *
 */

inline double Second1D::value() const {
  return v_;
}

inline void Second1D::setValue(double v) {
  v_ = v;
}

inline double Second1D::first() const {
  return d_;
}

inline void Second1D::setFirst(double d) {
  d_ = d;
}

inline double Second1D::second() const {
  return h_;
}

inline void Second1D::setSecond(double h) {
  h_ = h;
}

inline void Second1D::setZero() {
  *this = Second1D();
}

/*
 *
 * Implementation of all differentiation rules.
 *
 */

inline const Second1D& Second1D::operator+=(const Second1D& rhs) {
  v_ += rhs.v_;
  d_ += rhs.d_;
  h_ += rhs.h_;
  return *this;
}

inline const Second1D& Second1D::operator-=(const Second1D& rhs) {
  v_ -= rhs.v_;
  d_ -= rhs.d_;
  h_ -= rhs.h_;
  return *this;
}

inline const Second1D& Second1D::operator*=(const Second1D& rhs) {
  h_ = h_ * rhs.v_ + 2 * d_ * rhs.d_ + v_ * rhs.h_;
  d_ = v_ * rhs.d_ + rhs.v_ * d_;
  v_ *= rhs.v_;
  return *this;
}

inline const Second1D& Second1D::operator/=(const Second1D& rhs) {
  h_ = (rhs.v_ * rhs.v_ * h_ - rhs.v_ * (2 * d_ * rhs.d_ + v_ * rhs.h_) + 2 * v_ * rhs.d_ * rhs.d_) /
       (rhs.v_ * rhs.v_ * rhs.v_);
  d_ = d_ / rhs.v_ - v_ / (rhs.v_ * rhs.v_) * rhs.d_;
  v_ /= rhs.v_;
  return *this;
}

inline const Second1D& Second1D::operator*=(double rhs) {
  h_ *= rhs;
  v_ *= rhs;
  d_ *= rhs;
  return *this;
}

inline const Second1D& Second1D::operator/=(double rhs) {
  h_ /= rhs;
  v_ /= rhs;
  d_ /= rhs;
  return *this;
}

inline Second1D Second1D::operator+(const Second1D& der) const {
  return {v_ + der.v_, d_ + der.d_, h_ + der.h_};
}

inline Second1D Second1D::operator-() {
  return {-v_, -d_, -h_};
}

inline Second1D Second1D::operator+() {
  return {v_, d_, h_};
}

inline Second1D Second1D::operator-(const Second1D& rhs) const {
  return {v_ - rhs.v_, d_ - rhs.d_, h_ - rhs.h_};
}

inline Second1D Second1D::operator*(double f) const {
  return {v_ * f, d_ * f, h_ * f};
}

inline Second1D Second1D::operator*(const Second1D& rhs) const {
  return {v_ * rhs.v_, v_ * rhs.d_ + rhs.v_ * d_, h_ * rhs.v_ + 2 * d_ * rhs.d_ + v_ * rhs.h_};
}

inline Second1D Second1D::operator/(double f) const {
  return {v_ / f, d_ / f, h_ / f};
}

inline Second1D Second1D::operator/(const Second1D& rhs) const {
  return {v_ / rhs.v_, d_ / rhs.v_ - v_ / (rhs.v_ * rhs.v_) * rhs.d_,
          (rhs.v_ * rhs.v_ * h_ - rhs.v_ * (2 * d_ * rhs.d_ + v_ * rhs.h_) + 2 * v_ * rhs.d_ * rhs.d_) /
              (rhs.v_ * rhs.v_ * rhs.v_)};
}

inline Second1D operator*(double f, const Second1D& der) {
  return der.operator*(f);
}

inline Second1D operator+(double f, const Second1D& h) {
  return {f + h.value(), h.first(), h.second()};
}

inline Second1D operator+(const Second1D& h, double f) {
  return operator+(f, h);
}

inline Second1D operator-(double f, const Second1D& h) {
  return {f - h.value(), -h.first(), -h.second()};
}

inline Second1D operator-(const Second1D& h, double f) {
  return {h.value() - f, h.first(), h.second()};
}

inline Second1D operator/(double v, const Second1D& der) {
  double x = der.value();
  double x2 = x * x;
  return {v / der.value(), (-v / x2) * der.first(), v / x2 * (2 * der.first() * der.first() / x - der.second())};
}

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine

#endif // AUTOMATICDIFFERENTIATION_SECOND1D_H
