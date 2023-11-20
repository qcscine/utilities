/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef AUTOMATICDIFFERENTIATION_FIRST1D_H
#define AUTOMATICDIFFERENTIATION_FIRST1D_H

#include "FirstBase.h"

namespace Scine {
namespace Utils {
namespace AutomaticDifferentiation {

/**
 * @class First1D First1D.h
 * @brief Class representing values in one dimension and allowing for the automatic calculation of first derivatives.
 */
class First1D : public FirstBase<double, First1D> {
 public:
  /**
   * @brief Constructor with no arguments given.
   */
  First1D();
  /**
   * @brief Constructor with two arguments.
   * @param v value of the function.
   * @param d derivative of the function.
   */
  First1D(double v, double d);
  void setZero();

  /**
   *
   * Add function name in the singular form to avoid ambiguities.
   *
   */
  double derivative() const {
    return derivatives();
  };
};

// This constructor initiates the object with value and derivative to Zero.
inline First1D::First1D() : First1D(0, 0) {
}

// This constructor initiates the object with value and derivative to double v and double d.
inline First1D::First1D(double v, double d) : FirstBase(v, d) {
}

inline void First1D::setZero() {
  *this = First1D();
}

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine

#endif // AUTOMATICDIFFERENTIATION_FIRST1D_H
