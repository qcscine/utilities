/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef AUTOMATICDIFFERENTIATION_FIRSTND_H
#define AUTOMATICDIFFERENTIATION_FIRSTND_H

#include "FirstBase.h"
#include <Eigen/Core>

namespace Scine {
namespace Utils {
namespace AutomaticDifferentiation {
/**
 * @class FirstND FirstND.h
 * @brief Class representing values in N dimensions and allowing for the automatic calculation of first derivatives
 *        in those N dimensions.
 */
class FirstND : public FirstBase<Eigen::MatrixXd, FirstND> {
 public:
  /**
   * @brief Constructor that takes no arguments.
   */
  FirstND();
  /**
   * @brief Constructor that takes the value of the function as a double and the derivatives as an Eigen::MatrixXd
   */
  FirstND(double v, Eigen::MatrixXd d);
  void setZero();

  /**
   * @brief Getter for the number of dimensions of the function.
   * @return int Number of dimensions.
   */
  int dimensions() const;
  /**
   * @brief Getter for the derivative.
   * @param int Index to the derivative that should be returned in the N dimensional derivatives object.
   * @return double
   */
  double derivative(int index) const;
};

// Constructor that takes no arguments and initializes the value and derivatives to Zero.
inline FirstND::FirstND() : FirstND(0, Eigen::MatrixXd(0, 0)) {
}

// Constructor that takes the value of the function as a double and the derivative as an Eigen::MatrixXd.
inline FirstND::FirstND(double v, Eigen::MatrixXd d) : FirstBase(v, std::move(d)) {
}

// Returns the size of the derivatives.
inline int FirstND::dimensions() const {
  return static_cast<int>(d_.size());
}

// Returns the derivative at dimension index.
inline double FirstND::derivative(int index) const {
  assert(0 <= index && index < d_.size());
  return d_(index);
}

inline void FirstND::setZero() {
  *this = FirstND();
}

} // namespace AutomaticDifferentiation
} // namespace Utils
} // namespace Scine

#endif // AUTOMATICDIFFERENTIATION_FIRSTND_H
