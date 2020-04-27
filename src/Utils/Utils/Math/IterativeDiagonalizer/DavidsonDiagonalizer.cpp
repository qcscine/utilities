/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "DavidsonDiagonalizer.h"
#include "IndirectPreconditionerEvaluator.h"
#include "IndirectSigmaVectorEvaluator.h"
#include <Utils/IO/Logger.h>
#include <Eigen/Eigenvalues>

namespace Scine {
namespace Utils {

template<class MatrixType, DavidsonBalancedType s>
DavidsonDiagonalizer<MatrixType, s>::DavidsonDiagonalizer(int eigenvaluesToCompute, int numberGuessVectors, int totalDimension)
  : eigenvaluesToCompute_(eigenvaluesToCompute), blockSize_(numberGuessVectors), subspaceDimension_(numberGuessVectors) {
  assert(eigenvaluesToCompute > 0 &&
         "Unintended behaviour: calculate negative amount of eigenvalues in Davidson diagonalizer.");
  if (eigenvaluesToCompute_ > blockSize_)
    throw InvalidDavidsonInputException();
  initialize(totalDimension);
}

template<class MatrixType, DavidsonBalancedType s>
void DavidsonDiagonalizer<MatrixType, s>::initialize(int dimension) {
  if (eigenvaluesToCompute_ > dimension) {
    eigenvaluesToCompute_ = dimension;
    Utils::Log::warning() << "Too many eigenvalues requested (" << eigenvaluesToCompute_ << "), they were set to "
                          << dimension << ".";
  }
  if (subspaceDimension_ > dimension) {
    subspaceDimension_ = eigenvaluesToCompute_;
    blockSize_ = subspaceDimension_;
    Utils::Log::warning() << "Initial subspace dimension too big (" << subspaceDimension_ << "), it was reduced to "
                          << eigenvaluesToCompute_ << ".";
  }
  assert(dimension >= blockSize_ && "Matrix dimension cannot be smaller than number of guess vectors!");
  guessVectors_ = Eigen::MatrixXd::Zero(dimension, subspaceDimension_);
  if (s == DavidsonBalancedType::standard) {
    guessVectors_.block(0, 0, subspaceDimension_, subspaceDimension_) =
        Eigen::MatrixXd::Identity(subspaceDimension_, subspaceDimension_);
  }
  else {
    srand(seed_);
    guessVectors_.block(0, 0, subspaceDimension_, subspaceDimension_) =
        Eigen::MatrixXd::Random(subspaceDimension_, subspaceDimension_);
  }
  maxDimension_ = dimension;
  maxIterations_ = dimension;
  sigmaVectorEvaluator_.reset();
  preconditionerEvaluator_.reset();
}

template<class MatrixType, DavidsonBalancedType s>
void DavidsonDiagonalizer<MatrixType, s>::setMaxIterations(int maxIterations) {
  maxIterations_ = maxIterations;
}

template<class MatrixType, DavidsonBalancedType s>
void DavidsonDiagonalizer<MatrixType, s>::setMatrixToDiagonalize(const MatrixType& matrix) {
  if (type_ == DavidsonDirectType::direct) { // if a direct method is in use, no matrix is needed.
    return;
  }
  matrixToDiagonalize_ = matrix;
  originalDiagonal_ = matrixToDiagonalize_.diagonal();
  initialize(matrixToDiagonalize_.cols());
}

template<class MatrixType, DavidsonBalancedType s>
void DavidsonDiagonalizer<MatrixType, s>::setGuess(const Eigen::MatrixXd& guessVectors) {
  guessVectors_ = guessVectors.leftCols(subspaceDimension_);
}

template<class MatrixType, DavidsonBalancedType s>
void DavidsonDiagonalizer<MatrixType, s>::setDavidsonType(DavidsonDirectType type) {
  type_ = type;
}

template<class MatrixType, DavidsonBalancedType s>
void DavidsonDiagonalizer<MatrixType, s>::setSeed(unsigned seed) {
  seed_ = seed;
}

template<class MatrixType, DavidsonBalancedType s>
void DavidsonDiagonalizer<MatrixType, s>::setSigmaVectorEvaluator(std::unique_ptr<SigmaVectorEvaluator<MatrixType>>&& evaluator) {
  sigmaVectorEvaluator_ = std::move(evaluator);
}

template<class MatrixType, DavidsonBalancedType s>
void DavidsonDiagonalizer<MatrixType, s>::setPreconditionerEvaluator(std::unique_ptr<PreconditionerEvaluator>&& evaluator) {
  preconditionerEvaluator_ = std::move(evaluator);
}

template<class MatrixType, DavidsonBalancedType s>
EigenContainer DavidsonDiagonalizer<MatrixType, s>::solve() {
  checkEvaluators();

  // Set initial size of guess matrix
  subspaceDimension_ = blockSize_;
  for (int i = 0; i < maxIterations_; ++i) {
    performIteration();
    Log::info() << "Iteration number: " << i << "\tSubspace dimension: " << subspaceDimension_ << ".";
    if (converged_) { // if full matrix still calculate everything
      Log::info() << "Converged in " << i << " iterations, dimensionality " << subspaceDimension_ << ".";
      return eigenPairs_;
    }
  }
  throw DavidsonNotConvergedException();
}

template<class MatrixType, DavidsonBalancedType s>
void DavidsonDiagonalizer<MatrixType, s>::performIteration() {
  // Orthogonalize the projection matrix
  if (s == DavidsonBalancedType::standard) {
    orthogonalizeSubspace(guessVectors_);
    guessVectors_.leftCols(subspaceDimension_).colwise().normalize();
  }

  const Eigen::MatrixXd& projector = guessVectors_.leftCols(subspaceDimension_);
  if (s == DavidsonBalancedType::balanced) {
    for (int i = basisOverlap_.cols(); i < subspaceDimension_; ++i) {
      basisOverlap_.conservativeResize(subspaceDimension_, subspaceDimension_);
      basisOverlap_.col(i) = projector.transpose() * projector.col(i);
      Eigen::VectorXd symmetricRow = basisOverlap_.col(i);
      basisOverlap_.row(i) = symmetricRow;
    }
  }

  Eigen::MatrixXd sigmaMatrix = sigmaVectorEvaluator_->evaluateSigmaVector(projector);
  Eigen::MatrixXd projectedMatrix = projector.transpose() * sigmaMatrix;

  Eigen::MatrixXd subspaceEVectors;
  Eigen::VectorXd subspaceEValues;

  if (s == DavidsonBalancedType::standard) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> subspaceDiagonalizer(projectedMatrix);
    subspaceEVectors = subspaceDiagonalizer.eigenvectors();
    subspaceEValues = subspaceDiagonalizer.eigenvalues();
  }
  else if (s == DavidsonBalancedType::balanced) {
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> subspaceDiagonalizer(projectedMatrix, basisOverlap_);
    subspaceEVectors = subspaceDiagonalizer.eigenvectors();
    subspaceEValues = subspaceDiagonalizer.eigenvalues();
  }

  // Check convergence:
  // Get a matrix with the first *eigenvaluesToCompute_* residual vectors as column.
  // basically, it calculated r_k = H*Psi_{k, approx} - h_k*Psi_{k, approx}
  // = (sigmaMatrix - h_k * bMatrix) * v_k
  // H: matrix to diagonalize
  // Psi_{k, approx} k-th approximated Ritz eigenvectors for the subspace
  // h_k: k-th approximated Ritz eigenvalue
  // sigmaMatrix: matrix with as columns the sigma vectors
  // bMatrix: matrix with as column the tentative vectors
  // v_k: k-th subspace eigenvector
  Eigen::MatrixXd residualVectors(sigmaMatrix.rows(), eigenvaluesToCompute_);
  for (int col = 0; col < eigenvaluesToCompute_; ++col) {
    residualVectors.col(col) = (sigmaMatrix - subspaceEValues(col) * projector) * subspaceEVectors.col(col);
  }

  Eigen::VectorXd residualNorms = residualVectors.colwise().norm();
  converged_ = true;
  for (int i = 0; i < residualNorms.size(); ++i) {
    if (residualNorms(i) > eigenvalueTol)
      converged_ = false;
  }
  if (converged_) {
    eigenPairs_ = {subspaceEValues.head(eigenvaluesToCompute_), projector * subspaceEVectors.leftCols(eigenvaluesToCompute_)};
    return;
  }
  expandSubspace(subspaceEValues, residualVectors, projector, residualNorms);
}

template<class MatrixType, DavidsonBalancedType s>
void DavidsonDiagonalizer<MatrixType, s>::expandSubspace(const Eigen::VectorXd& subspaceEValues,
                                                         const Eigen::MatrixXd& residualVectors, const Eigen::MatrixXd& projector,
                                                         const Eigen::VectorXd& residualNorms) {
  // look at each residual vector
  for (int i = 0; i < eigenvaluesToCompute_; ++i) {
    if (residualNorms(i) > eigenvalueTol) {
      Eigen::VectorXd correctionVector =
          residualVectors.col(i).cwiseProduct(-preconditionerEvaluator_->evaluate(subspaceEValues(i)));
      if (s == DavidsonBalancedType::standard) {
        Eigen::VectorXd orthogonalizedCorrectionVector(guessVectors_.rows());
        orthogonalizedCorrectionVector = orthogonalizeToSubspace(correctionVector.normalized(), projector);

        if (orthogonalizedCorrectionVector.norm() / correctionVector.norm() > /*correctionTol*/ 1e-3) {
          if (subspaceDimension_ < maxDimension_) {
            ++subspaceDimension_;
            guessVectors_.conservativeResize(Eigen::NoChange, subspaceDimension_);
            guessVectors_.col(subspaceDimension_ - 1) = orthogonalizedCorrectionVector.normalized();
          }
        }
      }
      else {
        if (subspaceDimension_ < maxDimension_) {
          ++subspaceDimension_;
          guessVectors_.conservativeResize(Eigen::NoChange, subspaceDimension_);
          guessVectors_.col(subspaceDimension_ - 1) = correctionVector;
        }
      }
    }
  }
}

template<class MatrixType, DavidsonBalancedType s>
void DavidsonDiagonalizer<MatrixType, s>::checkEvaluators() {
  if (!sigmaVectorEvaluator_ || !preconditionerEvaluator_) {
    switch (type_) {
      case DavidsonDirectType::standard:
        setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<MatrixType>>(matrixToDiagonalize_));
        setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(originalDiagonal_));
        break;
      case DavidsonDirectType::direct:
        throw InvalidSigmaVectorEvaluator();
      default:
        throw InvalidDavidsonTypeException();
    }
  }
}

template<class MatrixType, DavidsonBalancedType s>
void DavidsonDiagonalizer<MatrixType, s>::orthogonalizeSubspace(Eigen::MatrixXd& trialSpace) {
  int spaceDimension = trialSpace.rows();
  // Calculate modified Gram-Schmidt QR decomposition
  Eigen::MatrixXd Q(spaceDimension, subspaceDimension_);
  Eigen::MatrixXd R(subspaceDimension_, subspaceDimension_);

  for (int i = 0; i < subspaceDimension_; ++i) {
    R(i, i) = trialSpace.col(i).norm();
    Q.col(i) = trialSpace.col(i) / R(i, i);
    for (int j = i + 1; j < subspaceDimension_; ++j) {
      R(i, j) = Q.col(i).transpose() * trialSpace.col(j);
      trialSpace.col(j) -= Q.col(i) * R(i, j);
    }
  }
}
template<class MatrixType, DavidsonBalancedType s>
Eigen::VectorXd DavidsonDiagonalizer<MatrixType, s>::orthogonalizeToSubspace(const Eigen::VectorXd& vector,
                                                                             const Eigen::MatrixXd& subspace) const {
  Eigen::VectorXd orthonormalizedVector = vector;
  for (int basis = 0; basis < subspace.cols(); ++basis) {
    orthonormalizedVector -= subspace.col(basis).dot(vector) * subspace.col(basis) / subspace.col(basis).norm();
  }
  return orthonormalizedVector;
}

template class DavidsonDiagonalizer<Eigen::MatrixXd>;
template class DavidsonDiagonalizer<Eigen::SparseMatrix<double>>;
template class DavidsonDiagonalizer<Eigen::MatrixXd, DavidsonBalancedType::balanced>;
template class DavidsonDiagonalizer<Eigen::SparseMatrix<double>, DavidsonBalancedType::balanced>;

} // namespace Utils
} // namespace Scine
