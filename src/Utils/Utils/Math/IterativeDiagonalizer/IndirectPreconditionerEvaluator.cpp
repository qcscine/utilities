/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "IndirectPreconditionerEvaluator.h"

namespace Scine {
namespace Utils {

IndirectPreconditionerEvaluator::IndirectPreconditionerEvaluator(const Eigen::VectorXd& diagonal)
  : diagonalOfMatrixToDiagonalize_(diagonal) {
}

Eigen::VectorXd IndirectPreconditionerEvaluator::evaluate(double eigenvalue) const {
  Eigen::VectorXd result = Eigen::inverse((diagonalOfMatrixToDiagonalize_.array() - eigenvalue));
  for (int i = 0; i < diagonalOfMatrixToDiagonalize_.size(); ++i) {
    if (std::isinf(result(i)))
      result(i) = 1.;
  }
  return result;
}

} // namespace Utils
} // namespace Scine
