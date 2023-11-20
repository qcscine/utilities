/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SubspaceCollapser.h"
#include "DavidsonDiagonalizer.h"
#include "SubspaceOrthogonalizer.h"

namespace Scine {
namespace Utils {

SubspaceCollapser::SubspaceCollapser() = default;

void SubspaceCollapser::setMaxSubspaceDimension(int maxSubspaceDimension) {
  maxSubspaceDimension_ = maxSubspaceDimension;
}
void SubspaceCollapser::setEigenvaluesToCompute(int eigenvaluesToCompute) {
  eigenvaluesToCompute_ = eigenvaluesToCompute;
}

bool SubspaceCollapser::collapseNeeded(const EigenContainer& eigenPairs, int subspaceDimension,
                                       std::vector<int> notConvergedRoots) {
  // Second condition is needed to prevent a subspace collapse from happening just after another one
  // having thus that oldSubspaceDimension_ > guessVectors_.cols()
  if (subspaceDimension >= maxSubspaceDimension_ && oldSubspaceDimension_ != 0) {
    notConvergedRoots_ = std::move(notConvergedRoots);
    newRitzEigenvectors_ = eigenPairs.eigenVectors.leftCols(eigenvaluesToCompute_);
    return true;
  }
  oldSubspaceDimension_ = subspaceDimension;
  oldRitzEigenvectors_ = eigenPairs.eigenVectors.leftCols(eigenvaluesToCompute_);
  return false;
}

Eigen::MatrixXd SubspaceCollapser::getCollapsedOrthogonalSubspace() {
  Eigen::MatrixXd collapsedSpace(newRitzEigenvectors_.rows(), eigenvaluesToCompute_ + notConvergedRoots_.size());
  collapsedSpace.leftCols(eigenvaluesToCompute_) = newRitzEigenvectors_.leftCols(eigenvaluesToCompute_);

  int index = 0;

  for (int root : notConvergedRoots_) {
    SubspaceOrthogonalizer::orthonormalizeToSubspace(oldRitzEigenvectors_.col(root),
                                                     collapsedSpace.leftCols(eigenvaluesToCompute_ + index),
                                                     collapsedSpace.col(eigenvaluesToCompute_ + index));
    ++index;
  }

  oldSubspaceDimension_ = 0;
  return collapsedSpace;
}

Eigen::MatrixXd SubspaceCollapser::getCollapsedNonOrthogonalSubspace() {
  Eigen::MatrixXd collapsedSpace(newRitzEigenvectors_.rows(), eigenvaluesToCompute_ + notConvergedRoots_.size());
  collapsedSpace.leftCols(eigenvaluesToCompute_) = newRitzEigenvectors_.leftCols(eigenvaluesToCompute_);
  SubspaceOrthogonalizer::qrOrthogonalize(collapsedSpace.leftCols(eigenvaluesToCompute_), eigenvaluesToCompute_);

  int index = 0;

  for (int root : notConvergedRoots_) {
    SubspaceOrthogonalizer::orthonormalizeToSubspace(oldRitzEigenvectors_.col(root),
                                                     collapsedSpace.leftCols(eigenvaluesToCompute_ + index),
                                                     collapsedSpace.col(eigenvaluesToCompute_ + index));
    ++index;
  }

  oldSubspaceDimension_ = 0;
  return collapsedSpace;
}

int SubspaceCollapser::calculateSubspaceCollapserIterations(int requiredRoots, int convergedRoots, int maxDimension) {
  int notConvergedRoots = requiredRoots - convergedRoots;
  int newDimension = 0;
  if (notConvergedRoots < 5) {
    newDimension = std::min(maxDimension, std::max(requiredRoots + 20 * notConvergedRoots, 80));
  }
  else if (notConvergedRoots < 30) {
    newDimension = std::min(maxDimension, std::max(requiredRoots + 6 * notConvergedRoots, 100));
  }
  else {
    newDimension = std::min(maxDimension, requiredRoots + 4 * notConvergedRoots);
  }
  return newDimension;
}
} // namespace Utils
} // namespace Scine
