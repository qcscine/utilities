/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_SIGMAVECTOREVALUATOR_H
#define UTILSOS_SIGMAVECTOREVALUATOR_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {

enum class SpinTransition { singlet, triplet };
/**
 * @class SigmaVectorEvaluator @file SigmaVectorEvaluator.h
 * @brief Interface for the sigma vector evaluators in the Davidson iterative diagonalizer.
 */
template<class MatrixType>
class SigmaVectorEvaluator {
 public:
  virtual ~SigmaVectorEvaluator() = default;
  virtual Eigen::MatrixXd evaluateSigmaVector(const Eigen::MatrixXd& guessVectors) const = 0;
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_SIGMAVECTOREVALUATOR_H
