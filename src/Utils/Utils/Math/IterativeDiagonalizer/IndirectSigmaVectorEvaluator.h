/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_INDIRECTSIGMAVECTOREVALUATOR_H
#define UTILSOS_INDIRECTSIGMAVECTOREVALUATOR_H

#include "SigmaVectorEvaluator.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace Scine {
namespace Utils {
/**
 * @class IndirectSigmaVectorEvaluator @file IndirectSigmaVectorEvaluator
 * @brief This sigma vector evaluator calculates the sigma vector as described in
 *        Davidson's iterative method paper.
 * E. R. Davidson,
 * The iterative calculation of a few of the lowest eigenvalues and corresponding eigenvectors
 * of large real-symmetric matrices, J. Comput. Phys, 1975, 17, 87-94,
 * https://doi.org/10.1016/0021-9991(75)90065-0
 */
template<class MatrixType>
class IndirectSigmaVectorEvaluator final : public SigmaVectorEvaluator<MatrixType> {
 public:
  explicit IndirectSigmaVectorEvaluator(const MatrixType& matrix) : fullMatrixToDiagonalize_(matrix) {
  }
  ~IndirectSigmaVectorEvaluator() final = default;

  Eigen::MatrixXd evaluateSigmaVector(const Eigen::MatrixXd& guessVectors) const final;

 private:
  MatrixType fullMatrixToDiagonalize_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_INDIRECTSIGMAVECTOREVALUATOR_H
