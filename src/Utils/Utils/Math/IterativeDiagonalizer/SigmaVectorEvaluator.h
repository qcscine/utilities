/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_SIGMAVECTOREVALUATOR_H
#define UTILSOS_SIGMAVECTOREVALUATOR_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @brief Denotes the spin symmetry in CIS and TD-DFT-like calculations on restricted references.
 */
enum class SpinTransition { Singlet, Triplet };
/**
 * @class SigmaVectorEvaluator @file SigmaVectorEvaluator.h
 * @brief Interface for the sigma vector evaluators in the Davidson iterative diagonalizer.
 */
class SigmaVectorEvaluator {
 public:
  virtual ~SigmaVectorEvaluator() = default;
  virtual const Eigen::MatrixXd& evaluate(const Eigen::MatrixXd& guessVectors) const = 0;

  /**
   * @brief Allows for internal handling of new subspace dimension.
   * @param newSubspaceDimension The new guess vectors number after subspace collapse.
   */
  virtual void collapsed(int newSubspaceDimension) = 0;
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_SIGMAVECTOREVALUATOR_H
