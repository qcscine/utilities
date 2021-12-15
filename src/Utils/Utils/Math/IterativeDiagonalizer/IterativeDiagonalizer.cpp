/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "IterativeDiagonalizer.h"
#include "DiagonalizerSettings.h"
#include "PreconditionerEvaluator.h"
#include "SigmaVectorEvaluator.h"
#include <Core/Log.h>
#include <chrono>
#include <iomanip>
#include <numeric>

namespace Scine {
namespace Utils {

IterativeDiagonalizer::IterativeDiagonalizer(int /*eigenvaluesToCompute*/, int /*totalDimension*/) {
}

IterativeDiagonalizer::~IterativeDiagonalizer() = default;

auto IterativeDiagonalizer::settings() const -> const Settings& {
  return *settings_;
}

auto IterativeDiagonalizer::settings() -> Settings& {
  return *settings_;
}

void IterativeDiagonalizer::applySettings() {
  settings_->check(maxDimension_);

  maxSubspaceExpansion_ = settings_->getInt(initialGuessDimensionOption);

  eigenvaluesToCompute_ = settings_->getInt(numberOfRootsOption);
  rootConverged_ = std::vector<bool>(eigenvaluesToCompute_, false);
}

void IterativeDiagonalizer::createGuess() {
  srand(settings_->getInt(seedOption));
  subspaceDimension_ = maxSubspaceExpansion_;
  guessVectors_ = Eigen::MatrixXd::Zero(maxDimension_, subspaceDimension_);
  guessVectors_.block(0, 0, subspaceDimension_, subspaceDimension_) =
      Eigen::MatrixXd::Identity(subspaceDimension_, subspaceDimension_) +
      0.01 * Eigen::MatrixXd::Random(subspaceDimension_, subspaceDimension_);
  guessVectors_.col(0) += 1e-5 * Eigen::VectorXd::Random(maxDimension_);
  if (guess_) {
    int maxCols = std::min(subspaceDimension_, int(guess_->cols()));
    guessVectors_.block(0, 0, guess_->rows(), maxCols) = guess_->leftCols(maxCols);
  }
}

void IterativeDiagonalizer::initialize(int dimension) {
  maxDimension_ = dimension;
  applySettings();
}

void IterativeDiagonalizer::setGuess(boost::optional<Eigen::MatrixXd> guessVectors) {
  if (guessVectors && guessVectors->rows() != maxDimension_) {
    throw std::runtime_error("Wrong row dimension in guess vector.");
  }
  guess_ = std::move(guessVectors);
}

void IterativeDiagonalizer::setSigmaVectorEvaluator(std::shared_ptr<SigmaVectorEvaluator> evaluator) {
  sigmaVectorEvaluator_ = std::move(evaluator);
}

void IterativeDiagonalizer::setPreconditionerEvaluator(std::shared_ptr<PreconditionerEvaluator> evaluator) {
  preconditionerEvaluator_ = std::move(evaluator);
}

const Utils::EigenContainer& IterativeDiagonalizer::getEigenPairs() const {
  return ritzEstimate_;
}

const SigmaVectorEvaluator& IterativeDiagonalizer::getSigmaVectorEvaluator() const {
  return *sigmaVectorEvaluator_;
}

SigmaVectorEvaluator& IterativeDiagonalizer::getSigmaVectorEvaluator() {
  return *sigmaVectorEvaluator_;
}

const EigenContainer& IterativeDiagonalizer::solve(Core::Log& log) {
  checkEvaluators();

  applySettings();
  createGuess();

  auto startDiagonalization = std::chrono::system_clock::now();

  printHeader(log);

  for (iteration_ = 0; iteration_ < settings_->getInt(SettingsNames::maxDavidsonIterations); ++iteration_) {
    auto start = std::chrono::system_clock::now();

    performIteration(log);
    auto stop = std::chrono::system_clock::now();
    iterationTime_ = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    printIteration(log);
    if (converged_) { // if full matrix still calculate everything
      auto endDiagonalization = std::chrono::system_clock::now();
      log.output << "CONVERGED!" << Core::Log::endl;
      auto totalTime =
          std::chrono::duration_cast<std::chrono::milliseconds>(endDiagonalization - startDiagonalization).count();
      log.output << "Time needed: " << totalTime << " ms.\n" << Core::Log::endl;
      return ritzEstimate_;
    }
  }
  throw DiagonalizerNotConvergedException();
}

void IterativeDiagonalizer::checkEvaluators() {
  if (!sigmaVectorEvaluator_ || !preconditionerEvaluator_) {
    throw InvalidSigmaVectorEvaluator();
  }
}

void IterativeDiagonalizer::printHeader(Core::Log& log) const {
  log.output << Core::Log::endl;
  log.output << std::setw(1) << "" << std::string(111, '=') << Core::Log::nl;
  log.output << std::setw(2) << "#" << std::setw(108) << "" << std::setw(2) << "#" << Core::Log::nl;
  log.output << std::setw(2) << "#" << std::setw(18) << "Iteration" << std::setw(18) << "Dimension" << std::setw(18)
             << "Max Residual" << std::setw(18) << "Min Space Norm" << std::setw(18) << "Roots Converged"
             << std::setw(18) << "Time [ms]" << std::setw(2) << "#" << Core::Log::nl;
  log.output << std::setw(2) << "#" << std::setw(108) << "" << std::setw(2) << "#" << Core::Log::nl;
  log.output << std::setw(1) << "" << std::string(111, '=') << Core::Log::endl;
}

void IterativeDiagonalizer::printIteration(Core::Log& log) const {
  int rootsConverged = std::count(rootConverged_.begin(), rootConverged_.end(), true);
  log.output << std::setw(2) << "" << std::setw(18) << iteration_ << std::setw(18) << subspaceDimension_ << std::setw(18)
             << residualNorms_.maxCoeff() << std::setw(18) << guessVectors_.colwise().norm().minCoeff() << std::setw(18)
             << rootsConverged << std::setw(18) << iterationTime_ << std::setw(2) << "" << Core::Log::endl;
}

} // namespace Utils
} // namespace Scine
