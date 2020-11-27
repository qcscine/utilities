/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_VECTORDERIVATIVES3D_H
#define UTILS_MATH_VECTORDERIVATIVES3D_H

#include "Second3D.h"

namespace Scine {
namespace Utils {
namespace AutomaticDifferentiation {

/**
 * @class VectorDerivatives3D VectorDerivatives3D.h
 * @brief Class to handle 3D derivatives of 3D vectors.
 */
class VectorDerivatives3D {
 public:
  using Second3D = AutomaticDifferentiation::Second3D;

  /**
   * @brief Default Constructor.
   */
  VectorDerivatives3D();
  /**
   * @brief Constructor from three Second3D objects.
   */
  VectorDerivatives3D(Second3D x, Second3D y, Second3D z);

  /**
   * @brief Accessor function for x direction.
   */
  const Second3D& x() const {
    return x_;
  }
  /**
   * @brief Accessor function for y direction.
   */
  const Second3D& y() const {
    return y_;
  }
  /**
   * @brief Accessor function for z direction.
   */
  const Second3D& z() const {
    return z_;
  }

  // Common objects:
  static VectorDerivatives3D spatialVectorHessian3D(const Eigen::Vector3d& v);
  static VectorDerivatives3D spatialVectorHessian3DWithInverseDerivative(const Eigen::Vector3d& v);

  // Arithmetic functions:
  VectorDerivatives3D operator+(const VectorDerivatives3D& rhs);
  VectorDerivatives3D operator-() const;
  VectorDerivatives3D operator*(double f) const;
  VectorDerivatives3D operator*(const Second3D& rhs) const;

  // Other Functions
  Second3D dot(const VectorDerivatives3D& rhs) const;
  Second3D dot(const Eigen::Vector3d& rhs) const;
  VectorDerivatives3D cross(const VectorDerivatives3D& rhs) const;
  VectorDerivatives3D cross(const Eigen::Vector3d& rhs) const;
  Second3D norm() const;

 private:
  Second3D x_, y_, z_;
};

inline VectorDerivatives3D::VectorDerivatives3D() : x_(Second3D()), y_(Second3D()), z_(Second3D()) {
}

inline VectorDerivatives3D::VectorDerivatives3D(Second3D x, Second3D y, Second3D z) : x_(x), y_(y), z_(z) {
}

inline VectorDerivatives3D VectorDerivatives3D::spatialVectorHessian3D(const Eigen::Vector3d& v) {
  return {Second3D(v.x(), 1, 0, 0), Second3D(v.y(), 0, 1, 0), Second3D(v.z(), 0, 0, 1)};
}

inline VectorDerivatives3D VectorDerivatives3D::spatialVectorHessian3DWithInverseDerivative(const Eigen::Vector3d& v) {
  return {Second3D(v.x(), -1, 0, 0), Second3D(v.y(), 0, -1, 0), Second3D(v.z(), 0, 0, -1)};
}

inline VectorDerivatives3D VectorDerivatives3D::operator+(const VectorDerivatives3D& rhs) {
  return {x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_};
}

inline VectorDerivatives3D VectorDerivatives3D::operator-() const {
  return {-x_, -y_, -z_};
}

inline VectorDerivatives3D VectorDerivatives3D::operator*(double f) const {
  return {x_ * f, y_ * f, z_ * f};
}

inline VectorDerivatives3D VectorDerivatives3D::operator*(const Second3D& rhs) const {
  return {x_ * rhs, y_ * rhs, z_ * rhs};
}

inline VectorDerivatives3D::Second3D VectorDerivatives3D::dot(const VectorDerivatives3D& rhs) const {
  return x_ * rhs.x_ + y_ * rhs.y_ + z_ * rhs.z_;
}

inline VectorDerivatives3D::Second3D VectorDerivatives3D::dot(const Eigen::Vector3d& rhs) const {
  return x_ * rhs.x() + y_ * rhs.y() + z_ * rhs.z();
}

inline VectorDerivatives3D VectorDerivatives3D::cross(const VectorDerivatives3D& rhs) const {
  return {y_ * rhs.z_ - z_ * rhs.y_, z_ * rhs.x_ - x_ * rhs.z_, x_ * rhs.y_ - y_ * rhs.x_};
}

inline VectorDerivatives3D VectorDerivatives3D::cross(const Eigen::Vector3d& rhs) const {
  return {y_ * rhs.z() - z_ * rhs.y(), z_ * rhs.x() - x_ * rhs.z(), x_ * rhs.y() - y_ * rhs.x()};
}

inline VectorDerivatives3D::Second3D VectorDerivatives3D::norm() const {
  Second3D root = sqrt(x_ * x_ + y_ * y_ + z_ * z_);
  return root;
}

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_VECTORDERIVATIVES3D_H
