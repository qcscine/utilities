/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_SUBSPACECOLLAPSER_H
#define UTILSOS_SUBSPACECOLLAPSER_H

#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {
struct EigenContainer;
/**
 * @class SubspaceCollapser @file SubspaceCollapser.h
 * @brief Class responsible for the subspace collapse procedure in the Davidson algorithm.
 * This class is a member of a Davidson algorithm, which it monitors for the need of
 * collapsing the guess vector subspace.
 * For more details, please consult:
 *
 * C. W. Murray, S. C. Racine, E. R. Davidson,
 * Improved Algorithms for the Lowest Few Eigenvalues and Associated Eigenvectors of Large Matrices,
 * J. Comp. Phys, 103, 1991, 382-389.
 *
 * C. David Sherrill, Henry F. Schaefer,
 * The Configuration Interaction Method: Advances in Highly Correlated Approaches,
 * Advances in Quantum Chemistry, 34, 1999, 143-269.
 */
class SubspaceCollapser {
 public:
  /**
   * @brief Constructor.
   */
  SubspaceCollapser();

  /**
   * @brief Sets the maximal dimension attainable before collapsing the space.
   */
  void setMaxSubspaceDimension(int maxSubspaceDimension);
  /**
   * @brief Sets the number of required eigenvalues.
   */
  void setEigenvaluesToCompute(int eigenvaluesToCompute);
  /**
   * @brief Checks whether the collapse is needed and fills the Ritz eigenvector estimates.
   * @param eigenVectors The subspace eigenvectors corresponding to the roots searched.
   * @param subspaceDimension The current guess vectors dimension.
   * @return Whether a subspace collapse is needed.
   */
  bool collapseNeeded(const EigenContainer& eigenPairs, int subspaceDimension, std::vector<int> notConvergedRoots);

  /**
   * @brief Returns the collapsed subspace.
   */
  Eigen::MatrixXd getCollapsedOrthogonalSubspace();
  Eigen::MatrixXd getCollapsedNonOrthogonalSubspace();

  /**
   * @brief Calculates a good dimension of the subspace with an empirical formula.
   */
  static int calculateSubspaceCollapserIterations(int requiredRoots, int convergedRoots, int maxDimension);

 private:
  Eigen::MatrixXd oldRitzEigenvectors_, newRitzEigenvectors_;
  std::vector<int> notConvergedRoots_;
  int oldSubspaceDimension_{};
  int eigenvaluesToCompute_{};
  int maxSubspaceDimension_ = 150;
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_SUBSPACECOLLAPSER_H
