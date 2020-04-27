/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ConvergenceChecker.h"
#include "ScfMethod.h"

namespace Scine {
namespace Utils {

ConvergenceChecker::ConvergenceChecker(const ScfMethod& method) : m_(method), converged_(false), threshold_(1e-5) {
}

void ConvergenceChecker::checkConvergence() {
  std::swap(oldMatrix_, newMatrix_);
  updateDensityMatrix();

  if (newMatrix_.size() != oldMatrix_.size()) {
    converged_ = false;
    return;
  }

  auto rmsd = calculateRMSD();

  converged_ = rmsd <= threshold_;
}

void ConvergenceChecker::updateDensityMatrix() {
  newMatrix_ = m_.getDensityMatrix().restrictedMatrix();
}

double ConvergenceChecker::calculateRMSD() {
  return (newMatrix_ - oldMatrix_).norm();
}
} // namespace Utils
} // namespace Scine
