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
#include <chrono>
#include <iomanip>

namespace Scine {
namespace Utils {

using namespace Utils::AutomaticDifferentiation;

ScfMethod::ScfMethod(bool unrestrictedCalculationPossible, Utils::DerivativeOrder maximalOrder, bool orthogonalBasisSet)
  : LcaoMethod(unrestrictedCalculationPossible, maximalOrder, orthogonalBasisSet),
    converged(false),
    performedIterations_(0),
    maxIterations(100),
    convergenceAccelerator_(*this) {
}

ScfMethod::~ScfMethod() = default;

void ScfMethod::initialize() {
  LcaoMethod::initialize();
  verifyPesValidity();
  reinitializeDensityMatrixGuess();
}

void ScfMethod::addModifier(std::shared_ptr<ScfModifier> modifier, int priority) {
  modifier->setMethod(this);
  modifier->initialize();
  auto p = std::find_if(modifiers.begin(), modifiers.end(),
                        [modifier](const ModifierContainer::value_type& m) { return m.second == modifier; });
  if (p == modifiers.end()) {
    if (priority < 0) {
      priority = 0;
    }
    if (priority > 10) {
      priority = 10;
    }
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
  verifyPesValidity();
  onConvergedCalculationStarts();
  performedIterations_ = 0;

  calculateDensityIndependentQuantities(d);
  for (auto& m : modifiers) {
    m.second->onOverlapCalculated();
  }
  printHeader(log);

  // Perform the first iteration
  performIteration(d);
  performedIterations_++;
  convergenceChecker_.update(*this);
  printIteration(log);

  // Loop until convergence
  converged = false;
  while (!convergenceChecker_.converged() && performedIterations_ < maxIterations) {
    performIteration(d);
    convergenceChecker_.update(*this);
    performedIterations_++;
    printIteration(log);
  }
  converged = convergenceChecker_.converged();

  // Get the gradient etc.
  finalizeCalculation(d);
  for (auto& m : modifiers) {
    m.second->onCalculationFinalized();
  }

  computeEnergyAndDerivatives(d);
  printFooter(log);
}

/*
 * Invariant part that is performed in every iteration,
 */
void ScfMethod::performIteration(Utils::Derivative d) {
  auto start = std::chrono::system_clock::now();
  for (auto& m : modifiers) {
    m.second->onIterationStart();
  }

  calculateDensityDependentQuantities(d);
  assembleFockMatrix();
  for (auto& m : modifiers) {
    m.second->onFockCalculated();
  }

  solveEigenValueProblem();
  for (auto& m : modifiers) {
    m.second->onGEPSolved();
  }

  calculateOccupation();
  calculateDensity();
  for (auto& m : modifiers) {
    m.second->onDensityCalculated();
  }
  // For convergence criterion
  electronicEnergy_ = electronicPart_->calculateElectronicEnergy();
  auto stop = std::chrono::system_clock::now();
  iterationTime_ = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
}

void ScfMethod::resetConvergenceCheck() {
  convergenceChecker_.set(convergenceChecker_.get());
  convergenceChecker_.update(*this);
}

void ScfMethod::evaluateDensity(Utils::Derivative derivativeOrder) {
  calculateDensityIndependentQuantities(derivativeOrder);
  calculateDensityDependentQuantities(derivativeOrder);
  finalizeCalculation(derivativeOrder);
  computeEnergyAndDerivatives(derivativeOrder);
  performedIterations_ = 1;
  converged = true;
}

void ScfMethod::onConvergedCalculationStarts() {
  if (densityMatrix_.numberElectrons() != nElectrons_) {
    reinitializeDensityMatrixGuess();
  }
  densityMatrix_.setUnrestricted(unrestrictedCalculationRunning_);
  electronicOccupationGenerator_->newScfCycleStarted();
}

void ScfMethod::solveEigenValueProblem() {
  if (basisSetIsOrthogonal()) {
    if (unrestrictedCalculationRunning_) {
      LcaoUtils::solveUnrestrictedEigenvalueProblem(fockMatrix_, eigenvectorMatrix_, singleParticleEnergies_);
    }
    else {
      LcaoUtils::solveRestrictedEigenvalueProblem(fockMatrix_, eigenvectorMatrix_, singleParticleEnergies_);
    }
  }
  else {
    if (unrestrictedCalculationRunning_) {
      LcaoUtils::solveUnrestrictedGeneralizedEigenvalueProblem(fockMatrix_, overlapMatrix_, eigenvectorMatrix_,
                                                               singleParticleEnergies_);
    }
    else {
      LcaoUtils::solveRestrictedGeneralizedEigenvalueProblem(fockMatrix_, overlapMatrix_, eigenvectorMatrix_,
                                                             singleParticleEnergies_);
    }
  }
}

void ScfMethod::reinitializeDensityMatrixGuess() {
  densityMatrix_ = densityMatrixGuess_->calculateGuess();
  if (unrestrictedCalculationRunning() && densityMatrix_.restricted()) {
    densityMatrix_.setAlphaAndBetaFromRestrictedDensity();
  }
}

void ScfMethod::calculateDensityDependentQuantities(Utils::Derivative d) {
  Utils::DerivativeOrder order = Utils::DerivativeOrder::Zero;
  if (d == Utils::Derivative::First) {
    order = Utils::DerivativeOrder::One;
  }
  if (d == Utils::Derivative::SecondAtomic || d == Utils::Derivative::SecondFull) {
    order = Utils::DerivativeOrder::Two;
  }

  electronicPart_->calculateDensityDependentPart(order);
}

void ScfMethod::finalizeCalculation(Utils::Derivative d) {
  Utils::DerivativeOrder order = Utils::DerivativeOrder::Zero;
  if (d == Utils::Derivative::First) {
    order = Utils::DerivativeOrder::One;
  }
  if (d == Utils::Derivative::SecondAtomic || d == Utils::Derivative::SecondFull) {
    order = Utils::DerivativeOrder::Two;
  }

  electronicPart_->finalize(order);

  assembleFockMatrix();
  solveEigenValueProblem();
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

void ScfMethod::printHeader(Core::Log& log) const {
  const auto names = convergenceChecker_.getNames();
  const auto extraSigns = static_cast<int>(names.size() * 25);
  log.output << Core::Log::endl;

  log.output << std::setw(1) << "" << std::string(68 + extraSigns, '=') << Core::Log::nl;
  log.output << std::right << std::setw(39 + extraSigns / 2) << "SCF Block" << Core::Log::endl;
  log.output << std::fixed << Core::Log::endl;

  log.output << std::setw(1) << "" << std::string(68 + extraSigns, '=') << Core::Log::nl;
  log.output << std::setw(2) << "#" << std::setw(65 + extraSigns) << "" << std::setw(2) << "#" << Core::Log::nl;
  log.output << std::setw(2) << "#" << std::setw(15) << "Iteration" << std::setw(25) << "Electronic Energy [Ha]";
  for (const auto& name : names) {
    log.output << std::setw(25) << name;
  }
  log.output << std::setw(25) << "Time [ms]" << std::setw(2) << "#" << Core::Log::nl;
  log.output << std::setw(2) << "#" << std::setw(65 + extraSigns) << "" << std::setw(2) << "#" << Core::Log::nl;
  log.output << std::setw(1) << "" << std::string(68 + extraSigns, '=') << Core::Log::endl;
}

void ScfMethod::printIteration(Core::Log& log) const {
  log.output << std::fixed << std::setprecision(10) << std::setw(2) << "" << std::setw(15) << performedIterations_
             << std::setw(25) << electronicEnergy_;
  for (const boost::optional<double>& current : convergenceChecker_.getCurrentValues()) {
    if (current) {
      log.output << std::setw(25) << *current;
    }
    else {
      log.output << std::setw(25) << "N/D";
    }
  }
  log.output << std::setw(25) << std::setprecision(5) << iterationTime_ << std::setw(2) << "" << Core::Log::endl;
}

void ScfMethod::printFooter(Core::Log& log) const {
  const auto names = convergenceChecker_.getNames();
  const auto extraSigns = static_cast<int>(names.size() * 25);

  log.output << std::setw(1) << "" << std::string(68 + extraSigns, '=') << Core::Log::nl;

  log.output << std::setprecision(10) << std::fixed << Core::Log::endl << Core::Log::endl;

  std::string negationBit = converged ? "" : "NOT ";
  log.output << std::right << std::setw(45) << negationBit + "CONVERGED AFTER " << performedIterations_ << " ITERATIONS"
             << Core::Log::endl;
  LcaoMethod::printFooter(log);
}
} // namespace Utils
} // namespace Scine
