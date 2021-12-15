/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_INDIRECTPRECONDITIONEREVALUATOR_H
#define UTILSOS_INDIRECTPRECONDITIONEREVALUATOR_H

#include "PreconditionerEvaluator.h"

namespace Scine {
namespace Utils {

/**
 * @class IndirectPreconditionerEvaluator @file IndirectPreconditionerEvaluator
 * @brief This preconditioner evaluator returns the diagonal preconditioner for the Davidson algorithm.
 */
class IndirectPreconditionerEvaluator final : public PreconditionerEvaluator {
 public:
  explicit IndirectPreconditionerEvaluator(const Eigen::VectorXd& diagonal);
  ~IndirectPreconditionerEvaluator() final = default;
  /**
   * @brief Evaluates the preconditioner vector.
   * The preconditioner elements are p_k = (H - h_k*I)^{-1},
   * where H is the diagonal of the full matrix.
   * @param eigenvalues The current guess for the eigenvalue h_k
   */
  Eigen::VectorXd evaluate(const Eigen::VectorXd& vectorToPrecondition, double eigenvalue) const final;

 private:
  Eigen::VectorXd diagonalOfMatrixToDiagonalize_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_INDIRECTPRECONDITIONEREVALUATOR_H
