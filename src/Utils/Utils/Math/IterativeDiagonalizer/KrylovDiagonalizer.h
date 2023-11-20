/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_DAVIDSONDIAGONALIZER_H
#define UTILS_DAVIDSONDIAGONALIZER_H

#include "IterativeDiagonalizer.h"

namespace Scine {
namespace Utils {

class SubspaceCollapser;

/**
 * @class DavidsonDiagonalizer @file DavidsonDiagonalizer.h
 * @brief Davidson-Liu iterative diagonalizer.
 * Two specialized variants of the algorithm are implemented:
 * - Classical algorithm as implemented by Davidson
 *   - Davidson, E. R., J. Comput. Phys., 1975, 17, 87-94.
 *   - Liu, B., The simultaneous expansion method for the iterative solution of large
 *     real-symmetric matrices, Numerical Algorithms in Chemistry: Algebraic Methods,
 *     Technical Report LBL-8158, pp 49-53.
 *   - R. M. Parrish, E. G. Hohenstein, T. J. Martinez, "Balancing" the Block
 *     Davidson-Liu Algorithm, J. Chem. Theory Comput., 2016, 12, 3003-3007.
 * - Balanced version as described by Furche et al.:
 *   - R. M. Parrish, E. G. Hohenstein, T. J. Martinez, "Balancing" the Block
 *     Davidson-Liu Algorithm, J. Chem. Theory Comput., 2016, 12, 3003-3007.
 *   - F. Furche, B. T. Kull, B. D. Nguyen, J. Kwon, Accelerating molecular property
 *     with nonorthonormal Krylov space methods, J. Chem. Phys, 2016, 144.
 *
 * While the "vanilla" implementation should work as an allrounder, the balanced
 * specialization is particularly effective if used in combination with some
 * density screening algorithm. For example, in direct CIS calculations the
 * Fock matrix elements to compute are conditioned on a Cauchy-Schwarz screening.
 * The non-normality condition of the balanced algorithm causes the norm of the residuals
 * to steadily decrease as the iterations go. This correspond to a steady increase
 * of the amount of elements of the Fock matrix that can be safely be ignored and
 * the consequent acceleration of the sigma vector construction.
 */
class KrylovDiagonalizer : public IterativeDiagonalizer {
 public:
  using IterativeDiagonalizer::settings_;
  /**
   * @brief Constructor.
   * Guess vectors are initialized as the identity matrix or a random diagonally dominant matrix.
   * If this is not suitable for the problem at hand, custom guess vectors can be given through the
   * DavidsonDiagonalizer::setGuess() method.
   */
  KrylovDiagonalizer(int eigenvaluesToCompute, int totalDimension);
  ~KrylovDiagonalizer() override;
  /**
   * @brief Method to apply the settings in the settings object.
   * This must be called to enact the changes in the settings.
   */
  void applySettings() override;

 protected:
  virtual void callCollapserImpl() = 0;
  virtual auto eigenDecomposition(const Eigen::MatrixXd& projectedMatrix) const -> EigenContainer = 0;
  virtual void onIterationStart();
  virtual void onSigmaMatrixEvaluation(const Eigen::MatrixXd& projector);
  virtual void filterCorrectionVectors(const Eigen::MatrixXd& projector, Eigen::MatrixXd& newGuessVectors) const;

  /// Generic methods
  void performIteration(Core::Log& log) final;
  void collapse(const Eigen::MatrixXd& projector, Core::Log& log);
  void calculateResiduals(const Eigen::MatrixXd& sigmaMatrix, const Eigen::MatrixXd& projector);
  void checkConvergence() override;
  /** Get a matrix with the first *eigenvaluesToCompute_* residual vectors as column.
   * basically, it calculates r_k = H*Psi_{k, approx} - h_k*Psi_{k, approx}
   * H: matrix to diagonalize
   * Psi_{k, approx} k-th approximated Ritz eigenvectors for the subspace,
   * which is equal to H*b*eigenVector_{k}.
   * h_k: k-th approximated Ritz eigenvalue
   */
  void expandSubspace(const Eigen::MatrixXd& projector);
  void addVectorsToGuessBasis(const Eigen::MatrixXd& newVectors);

  Eigen::MatrixXd residualVectors_;
  EigenContainer subspaceEigenpairs_;
  std::vector<int> notConvergedRoots_;
  std::unique_ptr<SubspaceCollapser> subspaceCollapser_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_DAVIDSONDIAGONALIZER_H
