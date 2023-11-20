/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "KrylovDiagonalizer.h"
#include "../MathUtils.h"
#include "DiagonalizerSettings.h"
#include "PreconditionerEvaluator.h"
#include "SigmaVectorEvaluator.h"
#include "SubspaceCollapser.h"
#include "SubspaceOrthogonalizer.h"
#include <Core/Log.h>
#include <Eigen/Eigenvalues>
#include <numeric>

namespace Scine {
namespace Utils {

KrylovDiagonalizer::KrylovDiagonalizer(int eigenvaluesToCompute, int totalDimension)
  : IterativeDiagonalizer(eigenvaluesToCompute, totalDimension) {
  settings_ = std::make_unique<KrylovSettings>(eigenvaluesToCompute, totalDimension);
  initialize(totalDimension);
  notConvergedRoots_.resize(eigenvaluesToCompute_);
  std::iota(notConvergedRoots_.begin(), notConvergedRoots_.end(), 0);
}

void KrylovDiagonalizer::applySettings() {
  IterativeDiagonalizer::applySettings();
  subspaceCollapser_ = std::make_unique<SubspaceCollapser>();
  subspaceCollapser_->setMaxSubspaceDimension(settings_->getInt(subspaceCollapseDimensionOption));
  subspaceCollapser_->setEigenvaluesToCompute(eigenvaluesToCompute_);
}

void KrylovDiagonalizer::performIteration(Core::Log& log) {
  // Handle for, for instance, orthogonalizing
  onIterationStart();

  const Eigen::MatrixXd& projector = guessVectors_.leftCols(subspaceDimension_);
  // Handle for, for instance, calculating basis overlap
  onSigmaMatrixEvaluation(projector);
  // Sigma vectors evaluation + Hamiltonian projection
  const Eigen::MatrixXd& sigmaMatrix = sigmaVectorEvaluator_->evaluate(projector);
  Eigen::MatrixXd projectedMatrix = (projector.transpose() * sigmaMatrix).selfadjointView<Eigen::Lower>();

  // Solve eigenproblem in subspace
  subspaceEigenpairs_ = eigenDecomposition(projectedMatrix);

  calculateResiduals(sigmaMatrix, projector);

  checkConvergence();

  if (converged_) {
    ritzEstimate_ = {ritzEstimate_.eigenValues.head(eigenvaluesToCompute_),
                     ritzEstimate_.eigenVectors.leftCols(eigenvaluesToCompute_)};
    return;
  }
  // Collapse subspace
  collapse(projector, log);
}

inline void KrylovDiagonalizer::calculateResiduals(const Eigen::MatrixXd& sigmaMatrix, const Eigen::MatrixXd& projector) {
  ritzEstimate_.eigenVectors = projector * subspaceEigenpairs_.eigenVectors.leftCols(eigenvaluesToCompute_);
  ritzEstimate_.eigenValues = subspaceEigenpairs_.eigenValues.head(eigenvaluesToCompute_);
  // Calculate residuals
  // TODO: Allow for expansion of more than just not converged roots (f.i. 2 new vectors per root)
  residualVectors_.resize(sigmaMatrix.rows(), notConvergedRoots_.size());
  int index = 0;
  for (int notConvergedRoot : notConvergedRoots_) {
    residualVectors_.col(index) = sigmaMatrix * subspaceEigenpairs_.eigenVectors.col(notConvergedRoot);
    residualVectors_.col(index).noalias() -=
        ritzEstimate_.eigenVectors.col(notConvergedRoot) * ritzEstimate_.eigenValues(notConvergedRoot);
    ++index;
  }
}

inline void KrylovDiagonalizer::checkConvergence() {
  residualNorms_ = residualVectors_.colwise().norm();
  numberConvergedRoots_ = 0;
  // Maps the index in the residual norms to the actual root number
  std::vector<int> oldNotConvergedRoots = notConvergedRoots_;
  notConvergedRoots_.clear();

  // Keep track of converged roots
  for (int i = 0; i < residualVectors_.cols(); ++i) {
    if (residualNorms_(i) < settings_->getDouble(residualNormToleranceOption)) {
      rootConverged_[oldNotConvergedRoots[i]] = true;
    }
    else {
      rootConverged_[oldNotConvergedRoots[i]] = false;
      notConvergedRoots_.push_back(oldNotConvergedRoots[i]);
    }
  }
  numberConvergedRoots_ = std::count(rootConverged_.begin(), rootConverged_.end(), true);
  converged_ = std::count(rootConverged_.begin(), rootConverged_.end(), false) == 0;

  assert(numberConvergedRoots_ <= eigenvaluesToCompute_);
}

void KrylovDiagonalizer::collapse(const Eigen::MatrixXd& projector, Core::Log& log) {
  subspaceCollapser_->setMaxSubspaceDimension(
      SubspaceCollapser::calculateSubspaceCollapserIterations(eigenvaluesToCompute_, numberConvergedRoots_, maxDimension_));
  if (!subspaceCollapser_->collapseNeeded(ritzEstimate_, subspaceDimension_, notConvergedRoots_)) {
    expandSubspace(projector);
  }
  else {
    callCollapserImpl();
    // If total space is small, one-shot this
    if (guessVectors_.cols() >= maxDimension_) {
      guessVectors_ = Eigen::MatrixXd::Identity(maxDimension_, maxDimension_);
      subspaceDimension_ = maxDimension_;
    }
    sigmaVectorEvaluator_->collapsed(subspaceDimension_);
    std::fill(rootConverged_.begin(), rootConverged_.end(), false);
    numberConvergedRoots_ = 0;
    converged_ = false;
    notConvergedRoots_.resize(eigenvaluesToCompute_);
    std::iota(notConvergedRoots_.begin(), notConvergedRoots_.end(), 0);
    log.output << Core::Log::nl << "Subspace collapsed. New dimension is " << subspaceDimension_ << Core::Log::nl
               << Core::Log::endl;

    printHeader(log);
  }
}

inline void KrylovDiagonalizer::expandSubspace(const Eigen::MatrixXd& projector) {
  // look at each residual vector
  Eigen::MatrixXd newGuessVectors(guessVectors_.rows(), notConvergedRoots_.size());
  int index = 0;
  for (int notConvergedRoot : notConvergedRoots_) {
    newGuessVectors.col(index) =
        preconditionerEvaluator_->evaluate(residualVectors_.col(index), ritzEstimate_.eigenValues(notConvergedRoot));
    ++index;
  }
  filterCorrectionVectors(projector, newGuessVectors);

  addVectorsToGuessBasis(newGuessVectors);
}

void KrylovDiagonalizer::filterCorrectionVectors(const Eigen::MatrixXd& /*projector*/,
                                                 Eigen::MatrixXd& /*newGuessVectors*/) const {
}

void KrylovDiagonalizer::addVectorsToGuessBasis(const Eigen::MatrixXd& newVectors) {
  int vectorsAdded = newVectors.cols();
  if (subspaceDimension_ + vectorsAdded < maxDimension_) {
    subspaceDimension_ += vectorsAdded;
    guessVectors_.conservativeResize(Eigen::NoChange, subspaceDimension_);
    guessVectors_.rightCols(vectorsAdded) = newVectors;
  }
  else {
    guessVectors_.conservativeResize(Eigen::NoChange, maxDimension_);
    guessVectors_.rightCols(maxDimension_ - subspaceDimension_) = newVectors.leftCols(maxDimension_ - subspaceDimension_);
    subspaceDimension_ = maxDimension_;
  }
}

void KrylovDiagonalizer::onIterationStart() {
}

void KrylovDiagonalizer::onSigmaMatrixEvaluation(const Eigen::MatrixXd& /*projector*/) {
}

KrylovDiagonalizer::~KrylovDiagonalizer() = default;

} // namespace Utils
} // namespace Scine
