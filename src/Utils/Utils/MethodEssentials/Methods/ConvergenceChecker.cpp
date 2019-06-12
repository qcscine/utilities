/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ConvergenceChecker.h"
#include "SCFMethod.h"

namespace Scine {
namespace Utils {

ConvergenceChecker::ConvergenceChecker(const SCFMethod& method) : m_(method), converged_(false), threshold_(1e-5) {
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
  double squareSum = 0;
  for (int i = 0; i < newMatrix_.size(); i++) {
    double change = newMatrix_(i) - oldMatrix_(i);
    squareSum += change * change;
  }
  return std::sqrt(squareSum / newMatrix_.size());
}
} // namespace Utils
} // namespace Scine
