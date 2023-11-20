/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "DavidsonDiagonalizer.h"
#include "../MathUtils.h"
#include "DiagonalizerSettings.h"
#include "PreconditionerEvaluator.h"
#include "SubspaceOrthogonalizer.h"
#include <Core/Log.h>
#include <Eigen/Eigenvalues>

namespace Scine {
namespace Utils {

NonOrthogonalDavidson::NonOrthogonalDavidson(int eigenvaluesToCompute, int totalDimension)
  : KrylovDiagonalizer(eigenvaluesToCompute, totalDimension) {
}

NonOrthogonalDavidson::~NonOrthogonalDavidson() = default;

EigenContainer NonOrthogonalDavidson::eigenDecomposition(const Eigen::MatrixXd& projectedMatrix) const {
  auto algorithm = settings_->getString(gepAlgorithmForBalancedMethodOption);
  if (algorithm == "standard")
    return MathUtils::stabilizedGeneralizedEigendecomposition<MathUtils::GepAlgorithm::Standard>(
        projectedMatrix, basisOverlap_.selfadjointView<Eigen::Upper>());
  else if (algorithm == "cholesky")
    return MathUtils::stabilizedGeneralizedEigendecomposition<MathUtils::GepAlgorithm::Cholesky>(
        projectedMatrix, basisOverlap_.selfadjointView<Eigen::Upper>());
  else if (algorithm == "simultaneous_diag")
    return MathUtils::stabilizedGeneralizedEigendecomposition<MathUtils::GepAlgorithm::SimultaneousDiagonalization>(
        projectedMatrix, basisOverlap_.selfadjointView<Eigen::Upper>());
  else
    throw InvalidDiagonalizerInputException("Algorithm " + algorithm + " not available for stable GEP solution.");
}

inline void NonOrthogonalDavidson::onSigmaMatrixEvaluation(const Eigen::MatrixXd& projector) {
  decltype(basisOverlap_)::Index newCols = subspaceDimension_ - basisOverlap_.cols();
  basisOverlap_.conservativeResize(subspaceDimension_, subspaceDimension_);
  Eigen::MatrixXd rightCols = projector.rightCols(newCols);
  basisOverlap_.rightCols(newCols) = projector.transpose() * rightCols;
}

inline void NonOrthogonalDavidson::filterCorrectionVectors(const Eigen::MatrixXd& projector,
                                                           Eigen::MatrixXd& newGuessVectors) const {
  double absNorm =
      std::abs(newGuessVectors.col(newGuessVectors.cols() - 1).transpose() * projector.col(projector.cols() - 1));
  if (absNorm > 1.0 - settings_->getDouble(correctionToleranceOption)) {
    newGuessVectors.col(newGuessVectors.cols() - 1) = residualVectors_.col(newGuessVectors.cols() - 1);
  }
}

void NonOrthogonalDavidson::callCollapserImpl() {
  guessVectors_ = subspaceCollapser_->getCollapsedNonOrthogonalSubspace();
  subspaceDimension_ = guessVectors_.cols();
  basisOverlap_.resize(0, 0);
}

OrthogonalDavidson::OrthogonalDavidson(int eigenvaluesToCompute, int totalDimension)
  : KrylovDiagonalizer(eigenvaluesToCompute, totalDimension) {
}

OrthogonalDavidson::~OrthogonalDavidson() = default;

inline void OrthogonalDavidson::filterCorrectionVectors(const Eigen::MatrixXd& projector, Eigen::MatrixXd& newGuessVectors) const {
  Eigen::MatrixXd orthogonalizedCorrectionVector = newGuessVectors;
  Eigen::VectorXd correctionVectorsNorms = newGuessVectors.colwise().norm();

  int dimension = 0;
  double maxNorm = 0;
  Eigen::VectorXd maxNormOrthogonalizedVector;
  for (int i = 0; i < orthogonalizedCorrectionVector.cols(); ++i) {
    Eigen::VectorXd orthogonalizedVector;
    SubspaceOrthogonalizer::orthogonalizeToSubspace(orthogonalizedCorrectionVector.col(i), projector, orthogonalizedVector);

    double norm = orthogonalizedVector.norm() / correctionVectorsNorms(i);
    if (norm > maxNorm) {
      std::swap(norm, maxNorm);
      maxNormOrthogonalizedVector = orthogonalizedVector;
    }

    if (norm / correctionVectorsNorms(i) > settings_->getDouble(correctionToleranceOption)) {
      newGuessVectors.col(dimension++) = orthogonalizedVector;
    }
  }
  // Retain one at least if everything is thrown away, the one that corrects the most
  if (dimension == 0) {
    newGuessVectors.col(0) = maxNormOrthogonalizedVector;
    dimension = 1;
  }
  newGuessVectors.conservativeResize(Eigen::NoChange, dimension);
}

inline void OrthogonalDavidson::onIterationStart() {
  KrylovDiagonalizer::onIterationStart();
  SubspaceOrthogonalizer::qrOrthogonalize(guessVectors_, subspaceDimension_);
  guessVectors_.leftCols(subspaceDimension_).colwise().normalize();
}

EigenContainer OrthogonalDavidson::eigenDecomposition(const Eigen::MatrixXd& projectedMatrix) const {
  EigenContainer result;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> subspaceDiagonalizer(projectedMatrix);
  result.eigenVectors = subspaceDiagonalizer.eigenvectors();
  result.eigenValues = subspaceDiagonalizer.eigenvalues();
  return result;
}

void OrthogonalDavidson::callCollapserImpl() {
  guessVectors_ = subspaceCollapser_->getCollapsedOrthogonalSubspace();
  subspaceDimension_ = guessVectors_.cols();
}
} // namespace Utils
} // namespace Scine
