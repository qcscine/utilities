/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LcaoMethod.h"
#include "OverlapCalculator.h"
#include "StructureDependentInitializer.h"
#include <Core/Log.h>
#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Utils/Scf/LcaoUtils/AufbauPrincipleOccupationGenerator.h>
#include <Utils/Scf/LcaoUtils/DensityMatrixGenerator.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupationGenerator.h>
#include <Utils/Scf/LcaoUtils/HomoLumoGapCalculator.h>
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>
#include <Utils/Scf/MethodInterfaces/AdditiveElectronicContribution.h>
#include <Utils/Scf/MethodInterfaces/ElectronicContributionCalculator.h>
#include <Utils/Scf/MethodInterfaces/RepulsionCalculator.h>
#include <Utils/Typenames.h>
#include <iomanip>

namespace Scine {
namespace Utils {

using namespace AutomaticDifferentiation;

LcaoMethod::LcaoMethod(bool unrestrictedCalculationPossible, DerivativeOrder maximalOrder, bool orthogonalBasisSet)
  : SinglePointMethod(maximalOrder),
    nAOs_(0),
    unrestrictedCalculationPossible_(unrestrictedCalculationPossible),
    basisSetIsOrthogonal_(orthogonalBasisSet) {
  setElectronicOccupationGenerator(std::make_unique<LcaoUtils::AufbauPrincipleOccupationGenerator>());
}

LcaoMethod::~LcaoMethod() = default;

const Eigen::MatrixXd& LcaoMethod::getOverlapMatrix() const {
  return overlapMatrix_;
}

void LcaoMethod::setOverlapMatrix(const Eigen::MatrixXd& S) {
  assert(S.rows() == nAOs_ && S.cols() == nAOs_ &&
         "Dimensions of given overlap matrix do not correspond to the number of atomic orbitals.");
  overlapMatrix_ = S;
}

void LcaoMethod::setFockMatrix(SpinAdaptedMatrix F) {
  // assert(F.rows() == nAOs_ && F.cols() == nAOs_ && "Dimensions of given Fock matrix do not correspond to the number
  // of atomic orbitals."); TODO
  fockMatrix_ = std::move(F);
}

void LcaoMethod::setDensityMatrix(DensityMatrix P) {
  // assert(P.rows() == nAOs_ && P.cols() == nAOs_ && "Dimensions of given density matrix do not correspond to the
  // number of atomic orbitals."); TODO
  densityMatrix_ = std::move(P);
  if (unrestrictedCalculationRunning() && !densityMatrix_.unrestricted()) {
    densityMatrix_.setUnrestricted(true);
  }
}

void LcaoMethod::setEnergyWeightedDensityMatrix(const Eigen::MatrixXd& W) {
  assert(W.rows() == nAOs_ && W.cols() == nAOs_ &&
         "Dimensions of given energy-weighted density matrix do not correspond to the number of atomic orbitals.");
  energyWeightedDensityMatrix_ = W;
}

void LcaoMethod::setMolecularOrbitals(MolecularOrbitals C) {
  // assert(C.rows() == nAOs_ && C.cols() == nAOs_ && "Dimensions of given eigenvector matrix do not correspond to the
  // number of atomic orbitals."); TODO
  eigenvectorMatrix_ = std::move(C);
}

void LcaoMethod::setUnrestrictedCalculation(bool b) {
  if (b && !unrestrictedCalculationPossible()) {
    throw std::runtime_error("Not possible to run unrestricted calculations with this method.");
  }

  if (unrestrictedCalculationRunning_ && !b) { // stop
    unrestrictedCalculationRunning_ = false;
  }
  else if (!unrestrictedCalculationRunning_ && b) { // start
    unrestrictedCalculationRunning_ = true;
    if (!densityMatrix_.unrestricted()) {
      densityMatrix_.setAlphaAndBetaFromRestrictedDensity();
    }
    if (!occupation_.isUnrestricted()) {
      occupation_.makeUnrestricted();
    }
    if (!eigenvectorMatrix_.isUnrestricted()) {
      eigenvectorMatrix_.makeUnrestricted();
    }
  }
}

void LcaoMethod::setSpinMultiplicity(int s) {
  assert(s > 0 && "The spin multiplicity must be larger than zero.");
  spinMultiplicity_ = s;
  if (spinMultiplicity_ != 1) {
    setUnrestrictedCalculation(true);
  }
}

double LcaoMethod::getHomoLumoGap() const {
  return LcaoUtils::HomoLumoGapCalculator::calculate(singleParticleEnergies_, occupation_);
}

void LcaoMethod::resizeLcaoMethodMatrices() {
  fockMatrix_.resize(nAOs_);
  overlapMatrix_.resize(nAOs_, nAOs_);
  densityMatrix_.resize(nAOs_);
  energyWeightedDensityMatrix_.resize(nAOs_, nAOs_);
  eigenvectorMatrix_.invalidate();
}

void LcaoMethod::setElectronicOccupationGenerator(std::unique_ptr<LcaoUtils::ElectronicOccupationGenerator>&& electronicOccupationSetter) {
  electronicOccupationGenerator_ = move(electronicOccupationSetter);
  electronicOccupationGenerator_->setMethod(this);
}

void LcaoMethod::calculateOccupation() {
  occupation_ = electronicOccupationGenerator_->generateOccupation();
}

void LcaoMethod::calculateDensity() {
  densityMatrix_ = LcaoUtils::DensityMatrixGenerator::generate(occupation_, eigenvectorMatrix_);
}

void LcaoMethod::calculateEnergyWeightedDensity() {
  energyWeightedDensityMatrix_ =
      LcaoUtils::DensityMatrixGenerator::generateEnergyWeighted(occupation_, eigenvectorMatrix_, singleParticleEnergies_);
}

void LcaoMethod::verifyPesValidity() const {
  verifyChargeValidity();
  verifyMultiplicityValidity();
  verifyUnrestrictedValidity();
}

void LcaoMethod::verifyChargeValidity() const {
  // Verify that the charge is not too positive (not taking away more electrons than available)
  if (molecularCharge_ > nElectronsForUnchargedSpecies_) {
    throw std::runtime_error("Not enough electrons to accomodate the molecular charge " + std::to_string(molecularCharge_) + ".");
  }
  else if (nElectronsForUnchargedSpecies_ - molecularCharge_ > 2 * nAOs_) {
    throw std::runtime_error("Not enough orbitals to accomodate the molecular charge " + std::to_string(molecularCharge_) + ".");
  }
}

void LcaoMethod::verifyMultiplicityValidity() const {
  if (spinMultiplicity_ > nElectrons_ + 1) {
    throw std::runtime_error("The chosen spin multiplicity (" + std::to_string(spinMultiplicity_) +
                             ") is too large (not enough electrons).");
  }
  if (spinMultiplicity_ > 2 * nAOs_ - nElectrons_ + 1) {
    throw std::runtime_error("The chosen spin multiplicity (" + std::to_string(spinMultiplicity_) +
                             ") is too large (not enough orbitals).");
  }
  // Check that number of electrons and spin multiplicity are compatible
  if ((spinMultiplicity_ + nElectrons_) % 2 == 0) {
    throw std::runtime_error("The chosen spin multiplicity (" + std::to_string(spinMultiplicity_) +
                             ") is not compatible with the molecular charge (" + std::to_string(molecularCharge_) + ").");
  }
}

void LcaoMethod::verifyUnrestrictedValidity() const {
  if (spinMultiplicity_ > 1 && !unrestrictedCalculationRunning_) {
    throw std::runtime_error("The chosen spin multiplicity (" + std::to_string(spinMultiplicity_) +
                             ") requires an unrestricted calculation.");
  }
}

void LcaoMethod::calculateBondOrderMatrix() {
  if (basisSetIsOrthogonal())
    LcaoUtils::calculateOrthonormalBondOrderMatrix(bondOrders_, densityMatrix_, aoIndexes_);
  else
    LcaoUtils::calculateBondOrderMatrix(bondOrders_, densityMatrix_, overlapMatrix_, aoIndexes_);
}

void LcaoMethod::calculateAtomicCharges() {
  if (basisSetIsOrthogonal())
    LcaoUtils::calculateOrthonormalAtomicCharges(atomicCharges_, coreCharges_, densityMatrix_, aoIndexes_);
  else
    LcaoUtils::calculateMullikenAtomicCharges(atomicCharges_, coreCharges_, densityMatrix_, overlapMatrix_, aoIndexes_);
}

void LcaoMethod::computeEnergyAndDerivatives(Utils::Derivative d) {
  electronicEnergy_ = electronicPart_->calculateElectronicEnergy();
  repulsionEnergy_ = rep_->getRepulsionEnergy();
  energy_ = electronicEnergy_ + repulsionEnergy_;

  if (d == Derivative::First)
    calculateDerivatives<Utils::Derivative::First>(gradients_);
  else if (d == Derivative::SecondAtomic)
    calculateDerivatives<Utils::Derivative::SecondAtomic>(secondDerivatives_);
  else if (d == Derivative::SecondFull)
    calculateDerivatives<Utils::Derivative::SecondFull>(fullSecondDerivatives_);
}

template<Utils::Derivative O>
void LcaoMethod::calculateDerivatives(DerivativeContainerType<O>& derivatives) {
  derivatives.setZero();

  rep_->addRepulsionDerivatives(derivatives);
  electronicPart_->addDerivatives(derivatives);
}

void LcaoMethod::initialize() {
  initializer_->initialize(elementTypes_);
  aoIndexes_ = initializer_->getAtomsOrbitalsIndexes();
  nAOs_ = aoIndexes_.getNAtomicOrbitals();
  nElectronsForUnchargedSpecies_ = initializer_->getNumberElectronsForUnchargedSpecies();
  nElectrons_ = nElectronsForUnchargedSpecies_ - molecularCharge_;
  coreCharges_ = initializer_->getCoreCharges();
  unrestrictedCalculationPossible_ = initializer_->unrestrictedCalculationPossible();

  overlapCalculator_->reset();
  electronicPart_->initialize();
  rep_->initialize();

  if (!unrestrictedCalculationPossible_)
    setUnrestrictedCalculation(false);

  resizeLcaoMethodMatrices();
  resizeRealTimeMethodMembers();
}

void LcaoMethod::calculateDensityIndependentQuantities(Utils::Derivative d) {
  DerivativeOrder order = DerivativeOrder::Zero;
  if (d == Derivative::First)
    order = DerivativeOrder::One;
  if (d == Derivative::SecondAtomic || d == Derivative::SecondFull)
    order = DerivativeOrder::Two;

  overlapCalculator_->calculateOverlap(order);
  overlapMatrix_ = overlapCalculator_->getOverlap().getMatrixXd();
  electronicPart_->calculateDensityIndependentPart(order);
  rep_->calculateRepulsion(order);
}

void LcaoMethod::assembleFockMatrix() {
  fockMatrix_ = electronicPart_->getMatrix();
}

void LcaoMethod::calculate(Utils::Derivative d, Core::Log& /*log*/) {
  verifyPesValidity();

  calculateDensityIndependentQuantities(d);
  assembleFockMatrix();
  LcaoUtils::solveRestrictedGeneralizedEigenvalueProblem(fockMatrix_, overlapMatrix_, eigenvectorMatrix_,
                                                         singleParticleEnergies_);

  calculateOccupationAndDensity();
  calculateBondOrderMatrix();
  calculateAtomicCharges();

  computeEnergyAndDerivatives(d);
}

void LcaoMethod::calculateOccupationAndDensity() {
  calculateOccupation();
  calculateDensity();
  if (!basisSetIsOrthogonal()) {
    calculateEnergyWeightedDensity();
  }
}

void LcaoMethod::addElectronicContribution(std::shared_ptr<AdditiveElectronicContribution> contribution) {
  if (contribution->isDensityDependent())
    electronicPart_->addDensityDependentElectronicContribution(std::move(contribution));
  else
    electronicPart_->addDensityIndependentElectronicContribution(std::move(contribution));
}

bool LcaoMethod::solvesOnlyOccupiedManifold() const {
  return solveOnlyOccupiedManifold_;
}

void LcaoMethod::setOnlyOccupiedManifoldToSolve(bool onlyOccupiedToSolve) {
  solveOnlyOccupiedManifold_ = onlyOccupiedToSolve;
}

void LcaoMethod::printFooter(Core::Log& log) const {
  log.output << std::setprecision(10) << std::fixed << Core::Log::endl << Core::Log::endl;
  log.output << std::setw(1) << "" << std::string(84, '=') << Core::Log::nl;
  log.output << std::setw(2) << "#" << std::setw(75) << "" << std::setw(8) << "#" << Core::Log::nl;
  log.output << std::setw(2) << "#" << std::setw(25) << "Electronic Energy";
  log.output << std::setw(25) << "Repulsion Energy" << std::setw(25) << "Total Energy" << std::setw(8) << "#"
             << Core::Log::nl;
  log.output << std::setw(2) << "#" << std::setw(22) << electronicEnergy_ << " Ha";
  log.output << std::setw(22) << repulsionEnergy_ << " Ha" << std::setw(22) << energy_ << " Ha" << std::setw(8) << "#"
             << Core::Log::nl;
  log.output << std::setw(2) << "#" << std::setw(75) << "" << std::setw(8) << "#" << Core::Log::nl;
  log.output << std::setw(1) << "" << std::string(84, '=') << Core::Log::endl;
  log.output << Core::Log::endl;
}

} // namespace Utils
} // namespace Scine
