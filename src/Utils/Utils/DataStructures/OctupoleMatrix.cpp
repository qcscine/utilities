/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "OctupoleMatrix.h"
#include "DipoleMatrix.h"

namespace Scine {
namespace Utils {

void OctupoleMatrix::reset(int dimension) {
  std::for_each(octupoleMatrices_.begin(), octupoleMatrices_.end(),
                [&](auto& octupoleMatrix) { octupoleMatrix.setDimension(dimension, dimension); });
}

Eigen::MatrixXd& OctupoleMatrix::operator[](int index) {
  return const_cast<Eigen::MatrixXd&>(const_cast<const OctupoleMatrix*>(this)->operator[](index));
}

const Eigen::MatrixXd& OctupoleMatrix::operator[](int index) const {
  return octupoleMatrices_.at(index).get0();
}

OctupoleMatrix::OctupoleMatrix(const OctupoleMatrix& rhs) = default;
OctupoleMatrix& OctupoleMatrix::operator=(const OctupoleMatrix& rhs) {
  if (this != &rhs) {
    octupoleMatrices_ = rhs.octupoleMatrices_;
  }
  return *this;
}
OctupoleMatrix::OctupoleMatrix() = default;
OctupoleMatrix::OctupoleMatrix(OctupoleMatrix&& rhs) noexcept = default;
OctupoleMatrix& OctupoleMatrix::operator=(OctupoleMatrix&& rhs) noexcept = default;
OctupoleMatrix::~OctupoleMatrix() = default;

} // namespace Utils
} // namespace Scine