/**
 * @file MathUtils.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INCLUDE_UTILS_MATH_UTILS_H
#define INCLUDE_UTILS_MATH_UTILS_H
#include "IterativeDiagonalizer/SpinAdaptedEigenContainer.h"
#include <Eigen/Eigenvalues>

namespace Scine {
namespace Utils {
namespace MathUtils {

/**
 * @brief Algorithms to use in the solution of the deflated generalized eigenvalue problem Ax=lBx.
 */
enum class GepAlgorithm {
  //! The standard Eigen::GeneralizedSelfAdjointEigenSolver. Uses a Cholesky LLT decomposition.
  Standard,
  //! An in-house implementation of the Cholesky LLT decomposition.
  Cholesky,
  //!* This algorithm works also for semidefinite B matrices.
  SimultaneousDiagonalization
};

struct StandardGepTag {};
struct CholeskyGepTag {};
struct SimultaneousDiagonalizationGepTag {};

template<GepAlgorithm gepAlgo = GepAlgorithm::Standard>
constexpr inline auto getTag() {
  return StandardGepTag{};
}
template<>
constexpr inline auto getTag<GepAlgorithm::Cholesky>() {
  return CholeskyGepTag{};
}
template<>
constexpr inline auto getTag<GepAlgorithm::SimultaneousDiagonalization>() {
  return SimultaneousDiagonalizationGepTag{};
}

/**
 * @brief Function to compute a generalized eigendecomposition.
 * Solves the generalized eigenvalue problem
 * Ax = lBx
 * using the standard Eigen::GeneralizedSelfAdjointEigenSolver.
 * the matrices A and B can also be deflated.
 * The eigenvalues are invariant under these transformations.
 *
 * @prec B is a positive definite matrix.
 * @prec A is a selfadjoint, real matrix.
 * @return The eigenvalues and eigenvectors, the solution of the problem.
 */
template<typename DerivedA, typename DerivedB>
inline auto eigendecompositionImpl(const DerivedA& transformedA, const DerivedB& preconditionedB, StandardGepTag)
    -> EigenContainer {
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(transformedA.template selfadjointView<Eigen::Upper>(),
                                                               preconditionedB.template selfadjointView<Eigen::Upper>());
  return EigenContainer{es.eigenvalues(), es.eigenvectors()};
}

/**
 * @brief Function to compute a stable generalized eigendecomposition.
 * Solves the generalized eigenvalue problem
 * Ax = lBx
 * in a stable way, using a preconditioning of the B matrix
 * B' = diag(B)^{-1/2} B diag(B)^{-1/2} ,
 * a Cholesky decomposition of the preconditioned B matrix
 * B' = LL^T,
 * and the definition of an auxiliary A matrix
 * A' = inv(L) * diag(B)^(-1/2) * A * diag(B)^(-1/2) * inv(L^T) ,
 * to solve the auxiliary Eigenvalue problem
 * A' * x' = O * x'
 * with
 * x' = L^T x .
 *
 * The original eigenvalues are recovered as
 * x = inv(L^T) x',
 * in which the inverse can be avoided by solution of the linear problem
 * L^T x = x' .
 * The eigenvalues are invariant under these transformations.
 *
 * @prec B is a positive definite matrix.
 * @prec A is a selfadjoint, real matrix.
 * @return The eigenvalues and eigenvectors, the solution of the problem.
 */
template<typename DerivedA, typename DerivedB>
inline auto eigendecompositionImpl(const DerivedA& A, const DerivedB& B, CholeskyGepTag) -> EigenContainer {
  // Compute the cholesky decomposition of matB = L L' = U'U
  Eigen::LLT<Eigen::MatrixXd> cholB(B.template selfadjointView<Eigen::Upper>());

  Eigen::MatrixXd matC = A.template selfadjointView<Eigen::Lower>();
  // compute C = inv(L) * d^(-1/2) * A * d^(-1/2) * inv(L')
  cholB.matrixL().solveInPlace<Eigen::OnTheLeft>(matC);
  cholB.matrixU().solveInPlace<Eigen::OnTheRight>(matC);
  // Solve C * x' = O * x'
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(matC);
  // transform back the eigen vectors x': x = d^(-1/2) * inv(U) * x'
  Eigen::MatrixXd eigenvectors = es.eigenvectors();
  cholB.matrixU().solveInPlace(eigenvectors);
  return EigenContainer{es.eigenvalues(), eigenvectors};
}

/**
 * @brief Solve the generalize eigenvalue problem via simultaneous diagonalization.
 * This method solves the generalized eigenvalue problem
 *   A * X = l * B * X
 * by finding the matrix X such that
 *   X^T * A * X = l
 * and
 *   X^T * B * X = I,
 * with X the eigenvectors of the matrix A and L a diagonal matrix having the eigenvalues as diagonal.
 * @param transformedA The A matrix, if B is preconditioned then it must be diag(B)^(-1/2) * A * diag(B)^(-1/2).
 * @param preconditionedB The B matrix, its eigenvalues can also be deflated by diag(B)^(-1/2) * B * diag(B)^(-1/2).
 * @return The eigenvalues/eigenvectors of the (possibly deflated) problem.
 */
template<typename DerivedA, typename DerivedB>
inline auto eigendecompositionImpl(const DerivedA& A, const DerivedB& B, SimultaneousDiagonalizationGepTag)
    -> EigenContainer {
  // Eigenvalue decomposition of the B matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B.template selfadjointView<Eigen::Lower>());
  // Count zero eigenvalues, calculate rank
  auto maxEigenvalue = es.eigenvalues().maxCoeff();
  auto threshold = maxEigenvalue / 1.0e8;
  int nZeroEigenvalues = 0;
  for (int i = 0; i < es.eigenvalues().size(); ++i)
    if (std::abs(es.eigenvalues()(i)) < threshold)
      ++nZeroEigenvalues;

  int rank = es.eigenvalues().size() - nZeroEigenvalues;
  // Take only eigenpairs for which the eigenvalue is != 0
  Eigen::VectorXd nonSingularEigenvalues = es.eigenvalues().tail(rank);
  Eigen::MatrixXd nonSingularEigenvectors = es.eigenvectors().rightCols(rank);

  // whiten the A matrix by A' = X'^T * A * X'
  // and X' is the whitening matrix X' = X_{B,r} * l_{B,r}^(-1/2)
  // X_{B,r}: eigenvectors of the B matrix, just the ones with non-zero eigenvalue
  // l_{B,r}^(-1/2): inverse sqrt of the non zero eigenvalues of B
  Eigen::MatrixXd whiteningTransformMatrix =
      nonSingularEigenvectors * nonSingularEigenvalues.array().sqrt().inverse().matrix().asDiagonal();
  Eigen::MatrixXd matC =
      whiteningTransformMatrix.transpose() * A.template selfadjointView<Eigen::Lower>() * whiteningTransformMatrix;

  // Solve the problem X''^T * A' * X'' = l
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esC(matC.selfadjointView<Eigen::Lower>());
  // Solution to the A * X = l * B * X, X = X' * X'' (rank(B) x rank(B) matrix)
  Eigen::MatrixXd transfFinal = whiteningTransformMatrix * esC.eigenvectors();

  return EigenContainer{esC.eigenvalues(), transfFinal};
}

/**
 * @brief Entry-point function to compute a generalized eigendecomposition with eigenvalues deflation.
 * Solves the generalized eigenvalue problem
 * Ax = lBx
 * in a stable way, using a preconditioning (deflation) of the B matrix
 * B' = diag(B)^{-1/2} B diag(B)^{-1/2}.
 * The matrix A must be transformed accordingly to
 * A' = diag(B)^{-1/2} B' diag(B)^{-1/2}.
 * The algorithm to perform the actual solution of the eigenproblem can be chosen to be
 * - Standard: the standard Eigen::GeneralizedSelfAdjointEigenSolver. Uses a Cholesky LLT decomposition.
 * - Cholesky: an in-house implementation of the Cholesky LLT decomposition.
 * - SimultaneousDiagonalization: works also for semidefinite B matrices.
 * After the solution of the deflated problem, the eigenvectors x are recovered as
 * x = d^(-1/2) * x',
 * with x' the solution of the deflated problem.
 */
template<GepAlgorithm gepAlgo = GepAlgorithm::Standard, typename DerivedA, typename DerivedB>
auto stabilizedGeneralizedEigendecomposition(const DerivedA& A, const DerivedB& B) -> EigenContainer {
  static_assert(gepAlgo == GepAlgorithm::Standard || gepAlgo == GepAlgorithm::Cholesky ||
                    gepAlgo == GepAlgorithm::SimultaneousDiagonalization,
                "Stable Generalize Eigenproblem algorithm not available.");
  // Precondition the basis overlap matrix to reduce the condition number.
  // S' = d^(-1/2) * S * d^(-1/2)
  Eigen::VectorXd diagB = B.diagonal().array().sqrt().inverse();
  Eigen::MatrixXd preconditionedB = diagB.asDiagonal() * Eigen::MatrixXd(B) * diagB.asDiagonal();
  Eigen::MatrixXd transformedA =
      (diagB.asDiagonal() * Eigen::MatrixXd(A) * diagB.asDiagonal()).template selfadjointView<Eigen::Lower>();

  EigenContainer result = eigendecompositionImpl(transformedA, preconditionedB, getTag<gepAlgo>());

  // Transform back the eigen vectors x': x = d^(-1/2) * x' due to B preconditioning
  result.eigenVectors.array().colwise() *= diagB.array();

  return result;
}

} // namespace MathUtils
} // namespace Utils
} // namespace Scine

#endif
