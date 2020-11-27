/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ScfMethod.h"
#include <Core/Log.h>
#include <Utils/Math/AutomaticDifferentiation/Second3D.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupationGenerator.h>
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>
#include <Utils/Scf/MethodInterfaces/DensityMatrixGuessCalculator.h>
#include <Utils/Scf/MethodInterfaces/ElectronicContributionCalculator.h>
#include <Utils/Scf/MethodInterfaces/OverlapCalculator.h>
#include <Utils/Scf/MethodInterfaces/RepulsionCalculator.h>
#include <Utils/Scf/MethodInterfaces/ScfModifier.h>

namespace Scine {
namespace Utils {

using namespace Utils::AutomaticDifferentiation;

ScfMethod::ScfMethod(bool unrestrictedCalculationPossible, Utils::DerivativeOrder maximalOrder, bool orthogonalBasisSet)
  : LcaoMethod(unrestrictedCalculationPossible, maximalOrder, orthogonalBasisSet),
    converged(false),
    performedIterations_(0),
    maxIterations(100),
    convergenceChecker_(*this),
    convergenceAccelerator_(*this) {
}

ScfMethod::~ScfMethod() = default;

void ScfMethod::initialize() {
  LcaoMethod::initialize();
  // Log is silent to prevent having to set it from outside and to avoid
  // unnecessary logging in the initialization phase.
  auto silentLog = Core::Log::silent();
  verifyPesValidity(silentLog);
  reinitializeDensityMatrixGuess();
}

void ScfMethod::addModifier(std::shared_ptr<ScfModifier> modifier, int priority) {
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

void ScfMethod::removeModifier(const std::shared_ptr<ScfModifier>& modifier) {
  for (auto it = modifiers.begin(); it != modifiers.end(); ++it) {
    if (it->second == modifier) {
      modifiers.erase(it);
      return;
    }
  }
}

void ScfMethod::calculate(Utils::Derivative d, Core::Log& log) {
  convergedCalculation(log, d);
}

/*
 * Perform a converged calculation
 */
void ScfMethod::convergedCalculation(Core::Log& log, Utils::Derivative d) {
  verifyPesValidity(log);
  onConvergedCalculationStarts();
  performedIterations_ = 0;

  calculateDensityIndependentQuantities(d);
  for (auto& m : modifiers)
    m.second->onOverlapCalculated();

  // Perform the first iteration
  performIteration(log, d);
  performedIterations_++;
  convergenceChecker_.checkConvergence();

  // Loop until convergence
  converged = false;
  while ((!convergenceChecker_.isConverged()) && performedIterations_ < maxIterations) {
    performIteration(log, d);
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
void ScfMethod::performIteration(Core::Log& log, Utils::Derivative d) {
  for (auto& m : modifiers)
    m.second->onIterationStart();

  calculateDensityDependentQuantities(d);
  assembleFockMatrix();
  for (auto& m : modifiers)
    m.second->onFockCalculated();

  solveEigenValueProblem(log);
  for (auto& m : modifiers)
    m.second->onGEPSolved();

  calculateOccupation();
  calculateDensity();
  for (auto& m : modifiers)
    m.second->onDensityCalculated();
}

void ScfMethod::resetConvergenceCheck() {
  convergenceChecker_.updateDensityMatrix();
}

void ScfMethod::evaluateDensity(Utils::Derivative derivativeOrder) {
  calculateDensityIndependentQuantities(derivativeOrder);
  calculateDensityDependentQuantities(derivativeOrder);
  assembleFockMatrix();
  finalizeCalculation(derivativeOrder);
  computeEnergyAndDerivatives(derivativeOrder);
  performedIterations_ = 1;
  converged = true;
}

void ScfMethod::onConvergedCalculationStarts() {
  if (densityMatrix_.numberElectrons() != nElectrons_) {
    reinitializeDensityMatrixGuess();
  }
  electronicOccupationGenerator_->newScfCycleStarted();
}

void ScfMethod::solveEigenValueProblem(Core::Log& log) {
  // Solve just for the occupied manifold, do not solve the virtual space.
  if (solvesOnlyOccupiedManifold()) {
    int alphaElectrons = 0, betaElectrons = 0;
    LcaoUtils::getNumberUnrestrictedElectrons(alphaElectrons, betaElectrons, nElectrons_, spinMultiplicity_);
    if (basisSetIsOrthogonal()) {
      if (unrestrictedCalculationRunning_)
        LcaoUtils::solveOccupiedUnrestrictedEigenvalueProblem(fockMatrix_, eigenvectorMatrix_, singleParticleEnergies_,
                                                              alphaElectrons, betaElectrons, log);
      else
        LcaoUtils::solveOccupiedRestrictedEigenvalueProblem(fockMatrix_, eigenvectorMatrix_, singleParticleEnergies_,
                                                            nElectrons_, log);
    }
    else {
      if (unrestrictedCalculationRunning_)
        LcaoUtils::solveOccupiedUnrestrictedGeneralizedEigenvalueProblem(
            fockMatrix_, overlapMatrix_, eigenvectorMatrix_, singleParticleEnergies_, alphaElectrons, betaElectrons, log);
      else
        LcaoUtils::solveOccupiedRestrictedGeneralizedEigenvalueProblem(fockMatrix_, overlapMatrix_, eigenvectorMatrix_,
                                                                       singleParticleEnergies_, nElectrons_, log);
    }
  }
  else { // Solve for the whole space.
    if (basisSetIsOrthogonal()) {
      if (unrestrictedCalculationRunning_)
        LcaoUtils::solveUnrestrictedEigenvalueProblem(fockMatrix_, eigenvectorMatrix_, singleParticleEnergies_);
      else
        LcaoUtils::solveRestrictedEigenvalueProblem(fockMatrix_, eigenvectorMatrix_, singleParticleEnergies_);
    }
    else {
      if (unrestrictedCalculationRunning_)
        LcaoUtils::solveUnrestrictedGeneralizedEigenvalueProblem(fockMatrix_, overlapMatrix_, eigenvectorMatrix_,
                                                                 singleParticleEnergies_);
      else
        LcaoUtils::solveRestrictedGeneralizedEigenvalueProblem(fockMatrix_, overlapMatrix_, eigenvectorMatrix_,
                                                               singleParticleEnergies_);
    }
  }
}

void ScfMethod::reinitializeDensityMatrixGuess() {
  densityMatrix_ = densityMatrixGuess_->calculateGuess();
  if (unrestrictedCalculationRunning() && densityMatrix_.restricted())
    densityMatrix_.setAlphaAndBetaFromRestrictedDensity();
}

void ScfMethod::calculateDensityDependentQuantities(Utils::Derivative d) {
  Utils::DerivativeOrder order = Utils::DerivativeOrder::Zero;
  if (d == Utils::Derivative::First)
    order = Utils::DerivativeOrder::One;
  if (d == Utils::Derivative::SecondAtomic || d == Utils::Derivative::SecondFull)
    order = Utils::DerivativeOrder::Two;

  electronicPart_->calculateDensityDependentPart(order);
}

void ScfMethod::finalizeCalculation(Utils::Derivative d) {
  Utils::DerivativeOrder order = Utils::DerivativeOrder::Zero;
  if (d == Utils::Derivative::First)
    order = Utils::DerivativeOrder::One;
  if (d == Utils::Derivative::SecondAtomic || d == Utils::Derivative::SecondFull)
    order = Utils::DerivativeOrder::Two;

  electronicPart_->finalize(order);

  calculateBondOrderMatrix();
  calculateAtomicCharges();
  if (!basisSetIsOrthogonal()) {
    calculateEnergyWeightedDensity();
  }
}

void ScfMethod::setScfMixer(scf_mixer_t mixer) {
  convergenceAccelerator_.setScfMixer(mixer);
}

scf_mixer_t ScfMethod::getScfMixer() const {
  return convergenceAccelerator_.getScfMixer();
}

DensityMatrix ScfMethod::getDensityMatrixGuess() const {
  return densityMatrixGuess_->calculateGuess();
}
} // namespace Utils
} // namespace Scine
