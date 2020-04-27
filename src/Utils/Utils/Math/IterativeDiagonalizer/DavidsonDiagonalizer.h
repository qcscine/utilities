/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_DAVIDSONDIAGONALIZER_H
#define UTILS_DAVIDSONDIAGONALIZER_H
#include "PreconditionerEvaluator.h"
#include "SigmaVectorEvaluator.h"
#include "SpinAdaptedEigenContainer.h"
#include <Eigen/Core>
#include <exception>
#include <map>
#include <memory>
namespace Scine {
namespace Utils {

class InvalidDavidsonInputException : public std::exception {
  const char* what() const noexcept final {
    return "Cannot calculate more eigenvalues than starting guess vectors!";
  }
};

class InvalidDavidsonTypeException : public std::exception {
  const char* what() const noexcept final {
    return "Only standard and direct types supported.";
  }
};

class InvalidSigmaVectorEvaluator : public std::exception {
  const char* what() const noexcept final {
    return "Empty sigma vector evaluator for Davidson iterative diagonalizer.";
  }
};

class DavidsonNotConvergedException : public std::exception {
  const char* what() const noexcept final {
    return "Davidson could not converge.";
  }
};
enum class DavidsonDirectType { standard, direct };
enum class DavidsonBalancedType { standard, balanced };

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
 * @tparam MatrixType The matrix type to use, can be either Eigen::MatrixXd or
 *         Eigen::SparseMatrix<double>.
 * @tparam s Determines if a standard (default) or balanced implementation is used.
 */
template<class MatrixType, DavidsonBalancedType s = DavidsonBalancedType::standard>
class DavidsonDiagonalizer {
 public:
  /**
   * @brief Constructor.
   * Guess vectors are initialized as the identity matrix. If this is not suitable
   * for the problem at hand, custom guess vectors can be given through the
   * DavidsonDiagonalizer::setGuess() method.
   */
  DavidsonDiagonalizer(int eigenvaluesToCompute, int numberGuessVectors, int totalDimension);
  /**
   * @brief Feeds in the matrix to be diagonalized.
   * @param matrix An Eigen::MatrixXd for which the first eigenvector/-value pairs are to be found.
   */
  void setMatrixToDiagonalize(const MatrixType& matrix);

  /**
   * @brief Sets the guess eigenvectors.
   * If this method is not called, the guess vectors are assumed to be the
   * identity matrix. This could not be ok for some problems (i.e. direct CIS method).
   * In these cases a custom set of guess eigenvectors must be given.
   */
  void setGuess(const Eigen::MatrixXd& guessVectors);
  /**
   * @brief Sets the sigma vector evaluator.
   * @param evaluator A class inheriting from Utils::SigmaVectorEvaluator.
   * For the indirect method an IndirectSigmaVectorEvaluator, which just calculates the matrix
   * multiplication of the guess vector and the underlying matrix.
   * For the direct method, a functor calculating the sigma vector without the need to store the whole
   * matrix.
   */
  void setSigmaVectorEvaluator(std::unique_ptr<SigmaVectorEvaluator<MatrixType>>&& evaluator);
  /**
   * @brief Sets the preconditioner evaluator.
   * @param evaluator A class inheriting from Utils::PreconditionerEvaluator.
   * For the indirect method an IndirectPreconditionerEvaluator, which just calculates the inverse
   * of the difference between the (approximate or exact) diagonal elements of the matrix to
   * diagonalize and the eigenvalues estimated so far.
   * For the direct method, a functor calculating the same quantity without the need to store the whole
   * matrix.
   */
  void setPreconditionerEvaluator(std::unique_ptr<PreconditionerEvaluator>&& evaluator);
  /**
   * @brief Set whether a standard or direct method has be used.
   * In that case no matrix must be given, just a functor for the direct calculation of the
   * sigma vectors.
   */
  void setDavidsonType(DavidsonDirectType type);
  /**
   * @brief Sets the maximal amount of iterations to reach convergence.
   */
  void setMaxIterations(int maxIterations);
  /**
   * @brief Sets the seed. Default value is 42.
   */
  void setSeed(unsigned seed);
  /**
   * @brief solves the eigenvalue problem iteratively for a subspace of the given matrix.
   * Depending on whether the standard or direct modality is chosen, the Davidson diagonalizer
   * uses the stored matrix or calculated on the fly the sigma matrix.
   */
  EigenContainer solve();

 protected:
  void initialize(int dimension);
  void checkEvaluators();
  void orthogonalizeSubspace(Eigen::MatrixXd& trialSpace);
  Eigen::VectorXd orthogonalizeToSubspace(const Eigen::VectorXd& vector, const Eigen::MatrixXd& subspace) const;
  void performIteration();
  /** Get a matrix with the first *eigenvaluesToCompute_* residual vectors as column.
   * basically, it calculates r_k = H*Psi_{k, approx} - h_k*Psi_{k, approx}
   * H: matrix to diagonalize
   * Psi_{k, approx} k-th approximated Ritz eigenvectors for the subspace,
   * which is equal to H*b*eigenVector_{k}.
   * h_k: k-th approximated Ritz eigenvalue
   */
  void expandSubspace(const Eigen::VectorXd& eigenvalues, const Eigen::MatrixXd& residualVectors,
                      const Eigen::MatrixXd& projector, const Eigen::VectorXd& residualNorms);

  MatrixType matrixToDiagonalize_;
  Eigen::MatrixXd guessVectors_;
  Eigen::MatrixXd basisOverlap_;
  Eigen::VectorXd originalDiagonal_;
  EigenContainer eigenPairs_;
  std::unique_ptr<SigmaVectorEvaluator<MatrixType>> sigmaVectorEvaluator_;
  std::unique_ptr<PreconditionerEvaluator> preconditionerEvaluator_;
  int eigenvaluesToCompute_{};
  int blockSize_{};
  int subspaceDimension_{};
  int maxDimension_{};
  int maxIterations_{};
  bool converged_{false};
  unsigned int seed_{42};
  DavidsonDirectType type_{DavidsonDirectType::standard};
  static constexpr const double eigenvalueTol = 7e-6;
  static constexpr const double correctionTol = 1e-5;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_DAVIDSONDIAGONALIZER_H
