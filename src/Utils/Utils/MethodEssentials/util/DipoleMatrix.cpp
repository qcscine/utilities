/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DipoleMatrix.h"

namespace Scine {
namespace Utils {

Eigen::MatrixXd& DipoleMatrix::x() {
  return xDipoleMatrix_;
}
Eigen::MatrixXd& DipoleMatrix::y() {
  return yDipoleMatrix_;
}
Eigen::MatrixXd& DipoleMatrix::z() {
  return zDipoleMatrix_;
}

Eigen::MatrixXd DipoleMatrix::x() const {
  return xDipoleMatrix_;
}

Eigen::MatrixXd DipoleMatrix::y() const {
  return yDipoleMatrix_;
}

Eigen::MatrixXd DipoleMatrix::z() const {
  return zDipoleMatrix_;
}

void DipoleMatrix::reset(int dimension) {
  xDipoleMatrix_ = yDipoleMatrix_ = zDipoleMatrix_ = Eigen::MatrixXd::Zero(dimension, dimension);
}

Eigen::MatrixXd& DipoleMatrix::operator[](int index) {
  if (index == 0)
    return x();
  else if (index == 1)
    return y();
  else if (index == 2)
    return z();
  throw InvalidDimensionSizeException();
}

Eigen::MatrixXd DipoleMatrix::operator[](int index) const {
  if (index == 0)
    return x();
  else if (index == 1)
    return y();
  else if (index == 2)
    return z();
  throw InvalidDimensionSizeException();
}

} // namespace Utils
} // namespace Scine
