/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "IndirectPreconditionerEvaluator.h"

namespace Scine {
namespace Utils {

IndirectPreconditionerEvaluator::IndirectPreconditionerEvaluator(const Eigen::VectorXd& diagonal)
  : diagonalOfMatrixToDiagonalize_(diagonal) {
}

Eigen::VectorXd IndirectPreconditionerEvaluator::evaluate(const Eigen::VectorXd& vectorToPrecondition, double eigenvalue) const {
  Eigen::VectorXd difference = eigenvalue - diagonalOfMatrixToDiagonalize_.array();
  Eigen::VectorXd preconditionedVector = vectorToPrecondition;
  for (int i = 0; i < diagonalOfMatrixToDiagonalize_.size(); ++i) {
    if (std::abs(difference(i)) < 1e-3) {
      preconditionedVector(i) = vectorToPrecondition(i);
    }
    else {
      preconditionedVector(i) /= difference(i);
    }
  }
  return preconditionedVector;
}

} // namespace Utils
} // namespace Scine
