/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "IndirectSigmaVectorEvaluator.h"
#include <exception>

namespace Scine {
namespace Utils {

template<class MatrixType>
const Eigen::MatrixXd& IndirectSigmaVectorEvaluator<MatrixType>::evaluate(const Eigen::MatrixXd& guessVectors) const {
  if (guessVectors.rows() != fullMatrixToDiagonalize_.cols())
    throw std::runtime_error("Dimensions of matrix to diagonalize and guess vector do not match.");
  cachedSigmaMatrix_ = fullMatrixToDiagonalize_.template selfadjointView<Eigen::Lower>() * guessVectors;
  return cachedSigmaMatrix_;
}

template<class MatrixType>
void IndirectSigmaVectorEvaluator<MatrixType>::collapsed(int /*newSubspaceDimension*/) {
}

template class IndirectSigmaVectorEvaluator<Eigen::MatrixXd>;
template class IndirectSigmaVectorEvaluator<Eigen::SparseMatrix<double>>;
} // namespace Utils
} // namespace Scine
