/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "QuadrupoleMatrix.h"
#include "DipoleMatrix.h"

namespace Scine {
namespace Utils {

void QuadrupoleMatrix::reset(int dimension) {
  std::for_each(quadrupoleMatrices_.begin(), quadrupoleMatrices_.end(),
                [&](auto& quadrupoleMatrix) { quadrupoleMatrix.setDimension(dimension, dimension); });
}

Eigen::MatrixXd& QuadrupoleMatrix::operator[](int index) {
  return const_cast<Eigen::MatrixXd&>(const_cast<const QuadrupoleMatrix*>(this)->operator[](index));
}

const Eigen::MatrixXd& QuadrupoleMatrix::operator[](int index) const {
  return quadrupoleMatrices_.at(index).get0();
}

QuadrupoleMatrix::QuadrupoleMatrix(const QuadrupoleMatrix& rhs) = default;
QuadrupoleMatrix& QuadrupoleMatrix::operator=(const QuadrupoleMatrix& rhs) {
  if (this != &rhs) {
    quadrupoleMatrices_ = rhs.quadrupoleMatrices_;
  }
  return *this;
}
QuadrupoleMatrix::QuadrupoleMatrix() = default;
QuadrupoleMatrix::QuadrupoleMatrix(QuadrupoleMatrix&& rhs) noexcept = default;
QuadrupoleMatrix& QuadrupoleMatrix::operator=(QuadrupoleMatrix&& rhs) noexcept = default;
QuadrupoleMatrix::~QuadrupoleMatrix() = default;
} // namespace Utils
} // namespace Scine