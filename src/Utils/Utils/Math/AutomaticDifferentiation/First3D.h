/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef AUTOMATICDIFFERENTIATION_FIRST3D_H
#define AUTOMATICDIFFERENTIATION_FIRST3D_H

#include "FirstBase.h"
#include <Eigen/Core>

namespace Scine {
namespace Utils {
namespace AutomaticDifferentiation {

/**
 * @class First3D First3D.h
 * @brief Class representing values in one dimensions and allowing for the automatic calculation of first derivatives.
 */
class First3D : public FirstBase<Eigen::Vector3d, First3D> {
 public:
  /**
   * @brief Constructor with no arguments given.
   */
  First3D();
  /**
   * @brief Constructor with the value and derivatives given as doubles.
   */
  First3D(double v, double d1, double d2, double d3);
  /**
   * @brief Constructor that takes the value as a double, but the derivative as an Eigen::Vector3d.
   */
  First3D(double v, Eigen::Vector3d d);
  void setZero();
};

// Constructor that takes no arguments and intializes the value and derivatives to Zero.
inline First3D::First3D() : First3D(0, Eigen::Vector3d::Zero()) {
}

// Constructor that takes the value and derivatives as doubles.
inline First3D::First3D(double v, double d1, double d2, double d3) : First3D(v, Eigen::Vector3d(d1, d2, d3)) {
}

// This constructor takes the derivatives as an Eigen::Vector3d object instead.
inline First3D::First3D(double v, Eigen::Vector3d d) : FirstBase(v, std::move(d)) {
}

inline void First3D::setZero() {
  *this = First3D();
}

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine

#endif // AUTOMATICDIFFERENTIATION_FIRST3D_H
