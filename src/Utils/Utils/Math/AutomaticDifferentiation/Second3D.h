/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef AUTOMATICDIFFERENTIATION_SECOND3D_H
#define AUTOMATICDIFFERENTIATION_SECOND3D_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {
namespace AutomaticDifferentiation {
/**
 * @class Second3D Second3D.h
 * @brief Class representing values in three dimensions and allowing for the automatic calculation of first and second
 *        derivatives.
 */
class Second3D {
 public:
  /**
   * @brief Constructor that takes no arguments and intializes value and derivatives to Zero.
   */
  Second3D();
  /**
   * @brief Constructor that takes the value and all derivatives as doubles.
   */
  Second3D(double v, double dx, double dy, double dz, double xx = 0, double yy = 0, double zz = 0, double xy = 0,
           double xz = 0, double yz = 0);

  /**
   * @brief Returns a Second3D object with first derivatives only multiplied by -1.
   * @return Second3D
   */
  Second3D opposite() const;

  /**
   * @brief Getter for the value.
   * @return double Value.
   */
  double value() const;
  /**
   * @brief Setter for the value.
   */
  void setValue(double v);
  /**
   * @brief Setter for first derivative as Eigen::Vector3d.
   */
  void setFirst3D(const Eigen::Ref<Eigen::Vector3d> d);
  /**
   * @brief Zero initializer.
   */
  void setZero();
  /**
   * @brief Getter for first derivative.
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d deriv() const;

  /**
   * @brief Setter for all the second derivatives.
   *
   *    XX refers to (dv)^2/(dx)^2, XY to (dv)^2/(dxdy), and so on...
   *
   */
  void setXX(double v);
  void setYY(double v);
  void setZZ(double v);
  void setXY(double v);
  void setYX(double v);
  void setXZ(double v);
  void setZX(double v);
  void setYZ(double v);
  void setZY(double v);

  /**
   * @brief Getter for dv/dx derivative.
   * @return double dv/dx
   */
  double dx() const;
  /**
   * @brief Getter for dv/dy derivative.
   * @return double dv/dy
   */
  double dy() const;
  /**
   * @brief Getter for dv/dz derivative.
   * @return double dv/dz
   */
  double dz() const;

  /**
   * @brief Getter for the second derivatives.
   * @return double
   */
  double XX() const;
  double YY() const;
  double ZZ() const;
  double XY() const;
  double YX() const;
  double XZ() const;
  double ZX() const;
  double YZ() const;
  double ZY() const;

  /**
   *
   * Definition of all differentiation rules.
   *
   */
  Second3D operator-() const;
  Second3D operator+() const;

  const Second3D& operator+=(const Second3D& rhs);
  const Second3D& operator-=(const Second3D& rhs);
  const Second3D& operator*=(const Second3D& rhs);
  const Second3D& operator/=(const Second3D& rhs);

  const Second3D& operator*=(double rhs);
  const Second3D& operator/=(double rhs);

  Second3D operator+(const Second3D& der) const;
  Second3D operator-(const Second3D& der) const;
  Second3D operator*(const Second3D& rhs) const;
  Second3D operator/(const Second3D& rhs) const;

  Second3D operator*(double f) const;
  Second3D operator/(double f) const;

 private:
  // value
  double v_;
  // first derivatives dv/dx, dv/dy and dv/dz
  double dx_, dy_, dz_;
  // second derivatives
  double xx_, yy_, zz_, xy_, xz_, yz_;
};

// Factor rule
Second3D operator*(double f, const Second3D& der);
// Deriving the square root function
Second3D sqrt(const Second3D& der);
// Deriving the exponential function
Second3D exp(const Second3D& der);
// Deriving the arccos function
Second3D arccos(const Second3D& der);

// Constructor that takes no arguments and intializes everything to Zero.
inline Second3D::Second3D() : Second3D(0, 0, 0, 0, 0, 0, 0, 0, 0, 0) {
}

// Constructor that takes everything as a double.
inline Second3D::Second3D(double v, double dx, double dy, double dz, double xx, double yy, double zz, double xy,
                          double xz, double yz)
  : v_(v), dx_(dx), dy_(dy), dz_(dz), xx_(xx), yy_(yy), zz_(zz), xy_(xy), xz_(xz), yz_(yz) {
}

/*
 *
 * Implementation of all differentiation rules
 *
 */

inline double Second3D::value() const {
  return v_;
}

inline void Second3D::setValue(double v) {
  v_ = v;
}

inline void Second3D::setFirst3D(const Eigen::Ref<Eigen::Vector3d> d) {
  dx_ = d.x();
  dy_ = d.y();
  dz_ = d.z();
}

inline void Second3D::setZero() {
  *this = Second3D();
}

inline Eigen::Vector3d Second3D::deriv() const {
  return {dx_, dy_, dz_};
}

inline void Second3D::setXX(double v) {
  xx_ = v;
}
inline void Second3D::setYY(double v) {
  yy_ = v;
}
inline void Second3D::setZZ(double v) {
  zz_ = v;
}
inline void Second3D::setXY(double v) {
  xy_ = v;
}
inline void Second3D::setYX(double v) {
  xy_ = v;
}
inline void Second3D::setXZ(double v) {
  xz_ = v;
}
inline void Second3D::setZX(double v) {
  xz_ = v;
}
inline void Second3D::setYZ(double v) {
  yz_ = v;
}
inline void Second3D::setZY(double v) {
  yz_ = v;
}

inline double Second3D::dx() const {
  return dx_;
}
inline double Second3D::dy() const {
  return dy_;
}
inline double Second3D::dz() const {
  return dz_;
}
inline double Second3D::XX() const {
  return xx_;
}
inline double Second3D::YY() const {
  return yy_;
}
inline double Second3D::ZZ() const {
  return zz_;
}
inline double Second3D::XY() const {
  return xy_;
}
inline double Second3D::YX() const {
  return xy_;
}
inline double Second3D::XZ() const {
  return xz_;
}
inline double Second3D::ZX() const {
  return xz_;
}
inline double Second3D::YZ() const {
  return yz_;
}
inline double Second3D::ZY() const {
  return yz_;
}

inline Second3D Second3D::opposite() const {
  return {v_, -dx_, -dy_, -dz_, xx_, yy_, zz_, xy_, xz_, yz_};
}

inline const Second3D& Second3D::operator+=(const Second3D& rhs) {
  v_ += rhs.v_;
  dx_ += rhs.dx_;
  dy_ += rhs.dy_;
  dz_ += rhs.dz_;
  xx_ += rhs.xx_;
  yy_ += rhs.yy_;
  zz_ += rhs.zz_;
  xy_ += rhs.xy_;
  xz_ += rhs.xz_;
  yz_ += rhs.yz_;
  return *this;
}

inline const Second3D& Second3D::operator-=(const Second3D& rhs) {
  v_ -= rhs.v_;
  dx_ -= rhs.dx_;
  dy_ -= rhs.dy_;
  dz_ -= rhs.dz_;
  xx_ -= rhs.xx_;
  yy_ -= rhs.yy_;
  zz_ -= rhs.zz_;
  xy_ -= rhs.xy_;
  xz_ -= rhs.xz_;
  yz_ -= rhs.yz_;
  return *this;
}

inline const Second3D& Second3D::operator*=(const Second3D& rhs) {
  xy_ = dx_ * rhs.dy_ + dy_ * rhs.dx_ + xy_ * rhs.v_ + v_ * rhs.xy_;
  xz_ = dx_ * rhs.dz_ + dz_ * rhs.dx_ + xz_ * rhs.v_ + v_ * rhs.xz_;
  yz_ = dy_ * rhs.dz_ + dz_ * rhs.dy_ + yz_ * rhs.v_ + v_ * rhs.yz_;
  xx_ = xx_ * rhs.v_ + 2 * dx_ * rhs.dx_ + v_ * rhs.xx_;
  yy_ = yy_ * rhs.v_ + 2 * dy_ * rhs.dy_ + v_ * rhs.yy_;
  zz_ = zz_ * rhs.v_ + 2 * dz_ * rhs.dz_ + v_ * rhs.zz_;
  dx_ = v_ * rhs.dx_ + rhs.v_ * dx_;
  dy_ = v_ * rhs.dy_ + rhs.v_ * dy_;
  dz_ = v_ * rhs.dz_ + rhs.v_ * dz_;
  v_ *= rhs.v_;
  return *this;
}

inline const Second3D& Second3D::operator/=(const Second3D& rhs) {
  double rhsV2 = rhs.v_ * rhs.v_;
  double invRhsV3 = 1.0 / (rhsV2 * rhs.v_);
  xx_ = (rhsV2 * xx_ - rhs.v_ * (2 * dx_ * rhs.dx_ + v_ * rhs.xx_) + 2 * v_ * rhs.dx_ * rhs.dx_) * invRhsV3;
  yy_ = (rhsV2 * yy_ - rhs.v_ * (2 * dy_ * rhs.dy_ + v_ * rhs.yy_) + 2 * v_ * rhs.dy_ * rhs.dy_) * invRhsV3;
  zz_ = (rhsV2 * zz_ - rhs.v_ * (2 * dz_ * rhs.dz_ + v_ * rhs.zz_) + 2 * v_ * rhs.dz_ * rhs.dz_) * invRhsV3;
  xy_ = (-rhs.v_ * (dx_ * rhs.dy_ + dy_ * rhs.dx_ + v_ * rhs.xy_) + xy_ * rhsV2 + 2 * v_ * rhs.dy_ * rhs.dx_) * invRhsV3;
  xz_ = (-rhs.v_ * (dx_ * rhs.dz_ + dz_ * rhs.dx_ + v_ * rhs.xz_) + xz_ * rhsV2 + 2 * v_ * rhs.dz_ * rhs.dx_) * invRhsV3;
  yz_ = (-rhs.v_ * (dz_ * rhs.dy_ + dy_ * rhs.dz_ + v_ * rhs.yz_) + yz_ * rhsV2 + 2 * v_ * rhs.dy_ * rhs.dz_) * invRhsV3;
  dx_ = dx_ / rhs.v_ - v_ / rhsV2 * rhs.dx_;
  dy_ = dy_ / rhs.v_ - v_ / rhsV2 * rhs.dy_;
  dz_ = dz_ / rhs.v_ - v_ / rhsV2 * rhs.dz_;
  v_ /= rhs.v_;
  return *this;
}

inline const Second3D& Second3D::operator*=(double rhs) {
  v_ *= rhs;
  dx_ *= rhs;
  dy_ *= rhs;
  dz_ *= rhs;
  xx_ *= rhs;
  yy_ *= rhs;
  zz_ *= rhs;
  xy_ *= rhs;
  xz_ *= rhs;
  yz_ *= rhs;
  return *this;
}

inline const Second3D& Second3D::operator/=(double rhs) {
  v_ /= rhs;
  dx_ /= rhs;
  dy_ /= rhs;
  dz_ /= rhs;
  xx_ /= rhs;
  yy_ /= rhs;
  zz_ /= rhs;
  xy_ /= rhs;
  xz_ /= rhs;
  yz_ /= rhs;
  return *this;
}

inline Second3D Second3D::operator+(const Second3D& der) const {
  return {v_ + der.v_,   dx_ + der.dx_, dy_ + der.dy_, dz_ + der.dz_, xx_ + der.xx_,
          yy_ + der.yy_, zz_ + der.zz_, xy_ + der.xy_, xz_ + der.xz_, yz_ + der.yz_};
}

inline Second3D Second3D::operator-() const {
  return {-v_, -dx_, -dy_, -dz_, -xx_, -yy_, -zz_, -xy_, -xz_, -yz_};
}

inline Second3D Second3D::operator+() const {
  return {v_, dx_, dy_, dz_, xx_, yy_, zz_, xy_, xz_, yz_};
}

inline Second3D Second3D::operator-(const Second3D& der) const {
  return {v_ - der.v_,   dx_ - der.dx_, dy_ - der.dy_, dz_ - der.dz_, xx_ - der.xx_,
          yy_ - der.yy_, zz_ - der.zz_, xy_ - der.xy_, xz_ - der.xz_, yz_ - der.yz_};
}

inline Second3D Second3D::operator*(const Second3D& rhs) const {
  return {v_ * rhs.v_,
          v_ * rhs.dx_ + rhs.v_ * dx_,
          v_ * rhs.dy_ + rhs.v_ * dy_,
          v_ * rhs.dz_ + rhs.v_ * dz_,
          xx_ * rhs.v_ + 2 * dx_ * rhs.dx_ + v_ * rhs.xx_,
          yy_ * rhs.v_ + 2 * dy_ * rhs.dy_ + v_ * rhs.yy_,
          zz_ * rhs.v_ + 2 * dz_ * rhs.dz_ + v_ * rhs.zz_,
          dx_ * rhs.dy_ + dy_ * rhs.dx_ + xy_ * rhs.v_ + v_ * rhs.xy_,
          dx_ * rhs.dz_ + dz_ * rhs.dx_ + xz_ * rhs.v_ + v_ * rhs.xz_,
          dy_ * rhs.dz_ + dz_ * rhs.dy_ + yz_ * rhs.v_ + v_ * rhs.yz_};
}

inline Second3D Second3D::operator*(double f) const {
  return {v_ * f, dx_ * f, dy_ * f, dz_ * f, xx_ * f, yy_ * f, zz_ * f, xy_ * f, xz_ * f, yz_ * f};
}

inline Second3D Second3D::operator/(double f) const {
  return {v_ / f, dx_ / f, dy_ / f, dz_ / f, xx_ / f, yy_ / f, zz_ / f, xy_ / f, xz_ / f, yz_ / f};
}

inline Second3D Second3D::operator/(const Second3D& rhs) const {
  double rhsV2 = rhs.v_ * rhs.v_;
  double invRhsV3 = 1.0 / (rhsV2 * rhs.v_);
  return {v_ / rhs.v_,
          dx_ / rhs.v_ - v_ / rhsV2 * rhs.dx_,
          dy_ / rhs.v_ - v_ / rhsV2 * rhs.dy_,
          dz_ / rhs.v_ - v_ / rhsV2 * rhs.dz_,
          (rhsV2 * xx_ - rhs.v_ * (2 * dx_ * rhs.dx_ + v_ * rhs.xx_) + 2 * v_ * rhs.dx_ * rhs.dx_) * invRhsV3,
          (rhsV2 * yy_ - rhs.v_ * (2 * dy_ * rhs.dy_ + v_ * rhs.yy_) + 2 * v_ * rhs.dy_ * rhs.dy_) * invRhsV3,
          (rhsV2 * zz_ - rhs.v_ * (2 * dz_ * rhs.dz_ + v_ * rhs.zz_) + 2 * v_ * rhs.dz_ * rhs.dz_) * invRhsV3,
          (-rhs.v_ * (dx_ * rhs.dy_ + dy_ * rhs.dx_ + v_ * rhs.xy_) + xy_ * rhsV2 + 2 * v_ * rhs.dy_ * rhs.dx_) * invRhsV3,
          (-rhs.v_ * (dx_ * rhs.dz_ + dz_ * rhs.dx_ + v_ * rhs.xz_) + xz_ * rhsV2 + 2 * v_ * rhs.dz_ * rhs.dx_) * invRhsV3,
          (-rhs.v_ * (dz_ * rhs.dy_ + dy_ * rhs.dz_ + v_ * rhs.yz_) + yz_ * rhsV2 + 2 * v_ * rhs.dy_ * rhs.dz_) * invRhsV3};
}

inline Second3D operator*(double f, const Second3D& der) {
  return der.operator*(f);
}

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine

#endif // AUTOMATICDIFFERENTIATION_SECOND3D_H
