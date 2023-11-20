/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DipoleMatrix.h"

namespace Scine {
namespace Utils {

void DipoleMatrix::reset(int dimension) {
  std::for_each(dipoleMatrices_.begin(), dipoleMatrices_.end(),
                [&](auto& dipoleMatrix) { dipoleMatrix.setDimension(dimension, dimension); });
}

Eigen::MatrixXd& DipoleMatrix::operator[](int index) {
  return const_cast<Eigen::MatrixXd&>(const_cast<const DipoleMatrix*>(this)->operator[](index));
}

const Eigen::MatrixXd& DipoleMatrix::operator[](int index) const {
  return dipoleMatrices_.at(index).get0();
}

DipoleMatrix::DipoleMatrix(const DipoleMatrix& rhs) = default;
DipoleMatrix& DipoleMatrix::operator=(const DipoleMatrix& rhs) {
  if (this != &rhs) {
    dipoleMatrices_ = rhs.dipoleMatrices_;
  }
  return *this;
}

DipoleMatrix::DipoleMatrix() = default;
DipoleMatrix::DipoleMatrix(DipoleMatrix&& rhs) noexcept = default;
DipoleMatrix& DipoleMatrix::operator=(DipoleMatrix&& rhs) noexcept = default;
DipoleMatrix::~DipoleMatrix() = default;

} // namespace Utils
} // namespace Scine
