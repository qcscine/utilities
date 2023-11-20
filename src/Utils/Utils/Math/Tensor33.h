/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_TENSOR33_H
#define UTILS_MATH_TENSOR33_H

#include <Eigen/Core>
#include <utility>

namespace Scine {
namespace Utils {

/**
 * @class Tensor33 Tensor33.h
 * @brief Simple algebra for a 3x3 tensor.
 */
class Tensor33 {
 public:
  /**
   * @brief Constructor from three Eigen::Vector3d vectors.
   */
  Tensor33(Eigen::Vector3d x, Eigen::Vector3d y, Eigen::Vector3d z);

  /**
   * @brief Accessor in x direction.
   */
  const Eigen::Vector3d& x() const {
    return x_;
  }
  /**
   * @brief Accessor in y direction.
   */
  const Eigen::Vector3d& y() const {
    return y_;
  }
  /**
   * @brief Accessor in z direction.
   */
  const Eigen::Vector3d& z() const {
    return z_;
  }

  // Mathematical operations
  Tensor33 operator+(const Tensor33& rhs) const;
  Tensor33 operator*(double f) const;

 private:
  Eigen::Vector3d x_, y_, z_;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
};

inline Tensor33::Tensor33(Eigen::Vector3d x, Eigen::Vector3d y, Eigen::Vector3d z)
  : x_(std::move(x)), y_(std::move(y)), z_(std::move(z)) {
}

inline Tensor33 Tensor33::operator+(const Tensor33& rhs) const {
  return {x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_};
}

inline Tensor33 Tensor33::operator*(double f) const {
  return {x_ * f, y_ * f, z_ * f};
}

inline Tensor33 tensor(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
  return {a.x() * b, a.y() * b, a.z() * b};
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_TENSOR33_H
