/**
 * @file IterativeDiagonalizer.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ITERATIVEDIAGONALIZER_H
#define UTILS_ITERATIVEDIAGONALIZER_H

#include "SpinAdaptedEigenContainer.h"
#include <Eigen/Core>
#include <boost/optional.hpp>
#include <exception>
#include <memory>

namespace Scine {
namespace Core {
struct Log;
} // namespace Core
namespace Utils {

class Settings;
class DiagonalizerSettings;
class SigmaVectorEvaluator;
class PreconditionerEvaluator;

/**
 * @brief Exception to throw if the signa vector evaluator is not valid, i.e. empty.
 */
class InvalidSigmaVectorEvaluator : public std::exception {
 public:
  const char* what() const noexcept final {
    return "Empty sigma vector evaluator or preconditioner for Davidson iterative diagonalizer.";
  }
};

/**
 * @brief Exception to throw if the diagonalizer does not converge.
 */
class DiagonalizerNotConvergedException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "Davidson could not converge.";
  }
};

/**
 * @brief Interface for iterative diagonalizers.
 * @class IterativeDiagonalizer
 */
class IterativeDiagonalizer {
 public:
  /**
   * @brief Constructor.
   * Guess vectors are initialized as the identity matrix or a random diagonally dominant matrix.
   * If this is not suitable for the problem at hand, custom guess vectors can be given through the
   * DavidsonDiagonalizer::setGuess() method.
   */
  IterativeDiagonalizer(int eigenvaluesToCompute, int totalDimension);
  virtual ~IterativeDiagonalizer();
  /**
   * @brief Const-getter for the Davidson settings.
   */
  auto settings() const -> const Settings&;
  /**
   * @brief Getter for the Davidson settings.
   */
  auto settings() -> Settings&;
  /**
   * @brief Method to apply the settings in the settings object.
   * This must be called to enact the changes in the settings.
   */
  virtual void applySettings();
  /**
   * @brief Sets the guess eigenvectors.
   * If this method is not called, the guess vectors are assumed to be the
   * identity matrix. This could not be ok for some problems (i.e. direct CIS method).
   * In these cases a custom set of guess eigenvectors must be given.
   */
  void setGuess(boost::optional<Eigen::MatrixXd> guessVectors);
  /**
   * @brief Sets the sigma vector evaluator.
   * @param evaluator A class inheriting from Utils::SigmaVectorEvaluator.
   * For the indirect method an IndirectSigmaVectorEvaluator, which just calculates the matrix
   * multiplication of the guess vector and the underlying matrix.
   * For the direct method, a functor calculating the sigma vector without the need to store the whole
   * matrix.
   */
  void setSigmaVectorEvaluator(std::shared_ptr<SigmaVectorEvaluator> evaluator);
  /**
   * @brief Sets the preconditioner evaluator.
   * @param evaluator A class inheriting from Utils::PreconditionerEvaluator.
   * For the indirect method an IndirectPreconditionerEvaluator, which just calculates the inverse
   * of the difference between the (approximate or exact) diagonal elements of the matrix to
   * diagonalize and the eigenvalues estimated so far.
   * For the direct method, a functor calculating the same quantity without the need to store the whole
   * matrix.
   */
  void setPreconditionerEvaluator(std::shared_ptr<PreconditionerEvaluator> evaluator);
  /**
   * @brief solves the eigenvalue problem iteratively for a subspace of the given matrix.
   * Depending on whether the standard or direct modality is chosen, the Davidson diagonalizer
   * uses the stored matrix or calculated on the fly the sigma matrix.
   */
  auto solve(Core::Log& log) -> const EigenContainer&;

  auto getEigenPairs() const -> const EigenContainer&;
  auto getSigmaVectorEvaluator() const -> const SigmaVectorEvaluator&;
  auto getSigmaVectorEvaluator() -> SigmaVectorEvaluator&;

 protected:
  void initialize(int dimension);
  void checkEvaluators();
  virtual void performIteration(Core::Log& log) = 0;
  virtual void checkConvergence() = 0;
  virtual void printHeader(Core::Log& log) const;
  virtual void printIteration(Core::Log& log) const;
  std::unique_ptr<DiagonalizerSettings> settings_;
  std::shared_ptr<SigmaVectorEvaluator> sigmaVectorEvaluator_;
  std::shared_ptr<PreconditionerEvaluator> preconditionerEvaluator_;
  boost::optional<Eigen::MatrixXd> guess_;
  Eigen::MatrixXd guessVectors_;
  EigenContainer ritzEstimate_;
  Eigen::VectorXd residualNorms_;
  double iterationTime_;
  std::vector<bool> rootConverged_;
  int numberConvergedRoots_{};
  int eigenvaluesToCompute_{};
  int subspaceDimension_{};
  int maxDimension_{};
  int maxSubspaceExpansion_{};
  int iteration_{};
  bool converged_{false};

 private:
  void createGuess();
};

} // namespace Utils
} // namespace Scine

#endif // Utils_ITERATIVEDIAGONALIZER_H
