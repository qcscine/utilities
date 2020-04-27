/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "IndirectSigmaVectorEvaluator.h"

namespace Scine {
namespace Utils {

template<class MatrixType>
Eigen::MatrixXd IndirectSigmaVectorEvaluator<MatrixType>::evaluateSigmaVector(const Eigen::MatrixXd& guessVectors) const {
  return fullMatrixToDiagonalize_.template selfadjointView<Eigen::Lower>() * guessVectors;
}

template class IndirectSigmaVectorEvaluator<Eigen::MatrixXd>;
template class IndirectSigmaVectorEvaluator<Eigen::SparseMatrix<double>>;
} // namespace Utils
} // namespace Scine
