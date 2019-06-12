/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SCFMethod.h"
#include <Utils/IO/Logger.h>
#include <Utils/Math/AutomaticDifferentiation/Second3D.h>
#include <Utils/MethodEssentials/Methods/DensityMatrixGuessCalculator.h>
#include <Utils/MethodEssentials/Methods/ElectronicContributionCalculator.h>
#include <Utils/MethodEssentials/Methods/OverlapCalculator.h>
#include <Utils/MethodEssentials/Methods/RepulsionCalculator.h>
#include <Utils/MethodEssentials/Methods/SCFModifier.h>
#include <Utils/MethodEssentials/util/LcaoUtil/ElectronicOccupationGenerator.h>
#include <Utils/MethodEssentials/util/LcaoUtil/LcaoUtil.h>

namespace Scine {
namespace Utils {

using namespace Utils::AutomaticDifferentiation;

SCFMethod::SCFMethod(bool unrestrictedCalculationPossible, Utils::derivOrder maximalOrder, bool orthogonalBasisSet)
  : LCAOMethod(unrestrictedCalculationPossible, maximalOrder, orthogonalBasisSet),
    converged(false),
    performedIterations_(0),
    maxIterations(100),
    convergenceChecker_(*this),
    convergenceAccelerator_(*this) {
}

SCFMethod::~SCFMethod() = default;

void SCFMethod::initialize() {
  LCAOMethod::initialize();
  reinitializeDensityMatrixGuess();
}

void SCFMethod::addModifier(std::shared_ptr<SCFModifier> modifier, int priority) {
  modifier->setMethod(this);
  modifier->initialize();
  auto p = std::find_if(modifiers.begin(), modifiers.end(),
                        [modifier](const ModifierContainer::value_type& m) { return m.second == modifier; });
  if (p == modifiers.end()) {
    if (priority < 0)
      priority = 0;
    if (priority > 10)
      priority = 10;
    modifiers.emplace(priority, modifier);
  }
}

void SCFMethod::removeModifier(const std::shared_ptr<SCFModifier>& modifier) {
  for (auto it = modifiers.begin(); it != modifiers.end(); ++it) {
    if (it->second == modifier) {
      modifiers.erase(it);
      return;
    }
  }
}

/*
 * Perform a converged calculation
 */
void SCFMethod::convergedCalculation(Utils::derivativeType d) {
  verifyPesValidity();
  onConvergedCalculationStarts();
  performedIterations_ = 0;

  calculateDensityIndependentQuantities(d);
  for (auto& m : modifiers)
    m.second->onOverlapCalculated();

  // Perform the first iteration
  performIteration(d);
  performedIterations_++;
  convergenceChecker_.checkConvergence();

  // Loop until convergence
  converged = false;
  while ((!convergenceChecker_.isConverged()) && performedIterations_ < maxIterations) {
    performIteration(d);
    convergenceChecker_.checkConvergence();
    performedIterations_++;
  }
  converged = convergenceChecker_.isConverged();

  // Get the gradient etc.
  finalizeCalculation(d);
  for (auto& m : modifiers)
    m.second->onCalculationFinalized();

  computeEnergyAndDerivatives(d);
}

/*
 * Invariant part that is performed in every iteration,
 */
void SCFMethod::performIteration(Utils::derivativeType d) {
  for (auto& m : modifiers)
    m.second->onIterationStart();

  calculateDensityDependentQuantities(d);
  assembleFockMatrix();
  for (auto& m : modifiers)
    m.second->onFockCalculated();

  solveEigenValueProblem();
  for (auto& m : modifiers)
    m.second->onGEPSolved();

  calculateOccupation();
  calculateDensity();
  for (auto& m : modifiers)
    m.second->onDensityCalculated();
}

void SCFMethod::resetConvergenceCheck() {
  convergenceChecker_.updateDensityMatrix();
}

void SCFMethod::evaluateDensity(Utils::derivativeType derivativeOrder) {
  calculateDensityIndependentQuantities(derivativeOrder);
  calculateDensityDependentQuantities(derivativeOrder);
  assembleFockMatrix();
  finalizeCalculation(derivativeOrder);
  computeEnergyAndDerivatives(derivativeOrder);
  performedIterations_ = 1;
  converged = true;
}

void SCFMethod::onConvergedCalculationStarts() {
  electronicOccupationGenerator_->newScfCycleStarted();
}
// TODO: Davidson for large systems
void SCFMethod::solveEigenValueProblem() {
  if (basisSetIsOrthogonal()) {
    if (unrestrictedCalculationRunning_)
      LcaoUtil::solveUnrestrictedEigenvalueProblem(fockMatrix_, eigenvectorMatrix_, singleParticleEnergies_);
    else
      LcaoUtil::solveRestrictedEigenvalueProblem(fockMatrix_, eigenvectorMatrix_, singleParticleEnergies_);
  }
  else {
    if (unrestrictedCalculationRunning_)
      LcaoUtil::solveUnrestrictedGeneralizedEigenvalueProblem(fockMatrix_, overlapMatrix_, eigenvectorMatrix_,
                                                              singleParticleEnergies_);
    else
      LcaoUtil::solveRestrictedGeneralizedEigenvalueProblem(fockMatrix_, overlapMatrix_, eigenvectorMatrix_,
                                                            singleParticleEnergies_);
  }
}

void SCFMethod::reinitializeDensityMatrixGuess() {
  densityMatrixGuess_->setNElectrons(nElectrons_);
  densityMatrix_ = densityMatrixGuess_->calculateGuess();
  if (unrestrictedCalculationRunning() && densityMatrix_.restricted())
    densityMatrix_.setAlphaAndBetaFromRestrictedDensity();
}

void SCFMethod::calculateDensityDependentQuantities(Utils::derivativeType d) {
  Utils::derivOrder order = Utils::derivOrder::zero;
  if (d == Utils::derivativeType::first)
    order = Utils::derivOrder::one;
  if (d == Utils::derivativeType::second_atomic || d == Utils::derivativeType::second_full)
    order = Utils::derivOrder::two;

  electronicPart_->calculateDensityDependentPart(order);
}

void SCFMethod::finalizeCalculation(Utils::derivativeType d) {
  Utils::derivOrder order = Utils::derivOrder::zero;
  if (d == Utils::derivativeType::first)
    order = Utils::derivOrder::one;
  if (d == Utils::derivativeType::second_atomic || d == Utils::derivativeType::second_full)
    order = Utils::derivOrder::two;

  electronicPart_->finalize(order);

  calculateBondOrderMatrix();
  calculateAtomicCharges();
  if (!basisSetIsOrthogonal()) {
    calculateEnergyWeightedDensity();
  }
}

void SCFMethod::setScfMixer(scf_mixer_t mixer) {
  convergenceAccelerator_.setScfMixer(mixer);
}

scf_mixer_t SCFMethod::getScfMixer() const {
  return convergenceAccelerator_.getScfMixer();
}

DensityMatrix SCFMethod::getDensityMatrixGuess() const {
  return densityMatrixGuess_->calculateGuess();
}
} // namespace Utils
} // namespace Scine
