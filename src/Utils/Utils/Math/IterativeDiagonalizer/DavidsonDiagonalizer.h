/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_NONORTHOGONALDAVIDSON_H
#define UTILS_NONORTHOGONALDAVIDSON_H

#include "KrylovDiagonalizer.h"

namespace Scine {
namespace Utils {

/**
 * @class NonOrthogonalDavidson @file DavidsonMethods.h
 * @brief Davidson-Liu iterative diagonalizer.
 * Balanced version as described by Furche et al.:
 *   - R. M. Parrish, E. G. Hohenstein, T. J. Martinez, "Balancing" the Block
 *     Davidson-Liu Algorithm, J. Chem. Theory Comput., 2016, 12, 3003-3007.
 *   - F. Furche, B. T. Kull, B. D. Nguyen, J. Kwon, Accelerating molecular property
 *     with nonorthonormal Krylov space methods, J. Chem. Phys, 2016, 144.
 *
 * The balanced specialization is particularly effective if used in combination with some
 * density screening algorithm. For example, in direct CIS calculations the
 * Fock matrix elements to compute are conditioned on a Cauchy-Schwarz screening.
 * The non-normality condition of the balanced algorithm causes the norm of the residuals
 * to steadily decrease as the iterations go. This correspond to a steady increase
 * of the amount of elements of the Fock matrix that can be safely be ignored and
 * the consequent acceleration of the sigma vector construction.
 */
class NonOrthogonalDavidson final : public KrylovDiagonalizer {
 public:
  /**
   * @brief Constructor.
   * Guess vectors are initialized as the identity matrix or a random diagonally dominant matrix.
   * If this is not suitable for the problem at hand, custom guess vectors can be given through the
   * DavidsonDiagonalizer::setGuess() method.
   */
  NonOrthogonalDavidson(int eigenvaluesToCompute, int totalDimension);
  ~NonOrthogonalDavidson() final;

 private:
  void onSigmaMatrixEvaluation(const Eigen::MatrixXd& projector) final;
  void callCollapserImpl() final;
  auto eigenDecomposition(const Eigen::MatrixXd& projectedMatrix) const -> EigenContainer final;
  /**
   * @brief Filters useless/redundant corrections in non-orthogonal Davidson
   * If the non-orthogonal davidson starts to give always the same correction
   * then don't use the preconditioned residual but just the residual.
   * So:
   * - Calulate the overlap of the new correction and the last vector in the space
   * as the problem presents itself exclusively for the convergence of the last root.
   * The rest is taken care with the subspace collapse.
   * - Look that it is greater than some threshold (I use 1 - the correction threshold
   * that is already a setting of the normal Davidson)
   * - If this is not the case, add the residual directly to the space instead
   * of the preconditioned residual
   */
  void filterCorrectionVectors(const Eigen::MatrixXd& projector, Eigen::MatrixXd& newGuessVectors) const final;
  Eigen::MatrixXd basisOverlap_;
};

/**
 * @class OrthogonalDavidson @file DavidsonMethods.h
 * @brief Davidson-Liu iterative diagonalizer.
 * Classical algorithm as implemented by Davidson
 *   - Davidson, E. R., J. Comput. Phys., 1975, 17, 87-94.
 *   - Liu, B., The simultaneous expansion method for the iterative solution of large
 *     real-symmetric matrices, Numerical Algorithms in Chemistry: Algebraic Methods,
 *     Technical Report LBL-8158, pp 49-53.
 *   - R. M. Parrish, E. G. Hohenstein, T. J. Martinez, "Balancing" the Block
 *     Davidson-Liu Algorithm, J. Chem. Theory Comput., 2016, 12, 3003-3007.
 *
 * "Vanilla" implementation, should work as an allrounder
 */
class OrthogonalDavidson final : public KrylovDiagonalizer {
 public:
  /**
   * @brief Constructor.
   * Guess vectors are initialized as the identity matrix or a random diagonally dominant matrix.
   * If this is not suitable for the problem at hand, custom guess vectors can be given through the
   * DavidsonDiagonalizer::setGuess() method.
   */
  OrthogonalDavidson(int eigenvaluesToCompute, int totalDimension);
  ~OrthogonalDavidson() final;

 private:
  void onIterationStart() final;
  void callCollapserImpl() final;
  auto eigenDecomposition(const Eigen::MatrixXd& projectedMatrix) const -> EigenContainer final;
  void filterCorrectionVectors(const Eigen::MatrixXd& projector, Eigen::MatrixXd& newGuessVectors) const final;
};
} // namespace Utils
} // namespace Scine

#endif // UTILS_DAVIDSONDIAGONALIZER_H
