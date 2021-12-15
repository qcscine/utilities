/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_PRECONDITIONEREVALUATOR_H
#define UTILSOS_PRECONDITIONEREVALUATOR_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @class PreconditionerEvaluator @file PreconditionerEvaluator.h
 * @brief Interface for the preconditioner evaluators in the Davidson iterative diagonalizer.
 */
class PreconditionerEvaluator {
 public:
  virtual ~PreconditionerEvaluator() = default;
  /**
   * @brief Evaluates the preconditioner vector.
   * The preconditioner elements are p_k = (H - h_k*I)^{-1},
   * where H is the exact or approximated matrix to diagonalize.
   * For diagonally dominant matrices H can be approximated with
   * D, the diagonal of H.
   */
  virtual Eigen::VectorXd evaluate(const Eigen::VectorXd& vectorToPrecondition, double eigenvalue) const = 0;
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_PRECONDITIONEREVALUATOR_H
