/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LcaoMethod.h"
#include "OverlapCalculator.h"
#include "StructureDependentInitializer.h"
#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <Utils/IO/Logger.h>
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

namespace Scine {
namespace Utils {

using namespace AutomaticDifferentiation;

LcaoMethod::LcaoMethod(bool unrestrictedCalculationPossible, derivOrder maximalOrder, bool orthogonalBasisSet)
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
  if (unrestrictedCalculationRunning() && !densityMatrix_.unrestricted())
    densityMatrix_.setUnrestricted(true);
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
    Log::error() << "Not possible to run unrestricted calculations with this method.";
    b = false;
  }

  if (unrestrictedCalculationRunning_ && !b) { // stop
    unrestrictedCalculationRunning_ = false;
  }
  else if (!unrestrictedCalculationRunning_ && b) { // start
    unrestrictedCalculationRunning_ = true;
    if (!densityMatrix_.unrestricted())
      densityMatrix_.setAlphaAndBetaFromRestrictedDensity();
    if (!occupation_.isUnrestricted())
      occupation_.makeUnrestricted();
    if (!eigenvectorMatrix_.isUnrestricted())
      eigenvectorMatrix_.makeUnrestricted();
  }
}

void LcaoMethod::setSpinMultiplicity(int s) {
  assert(s > 0 && "The spin multiplicity must be larger than zero.");
  spinMultiplicity_ = s;
}

double LcaoMethod::getHomoLumoGap() const {
  try {
    return LcaoUtils::HomoLumoGapCalculator::calculate(singleParticleEnergies_, occupation_);
  }
  catch (LcaoUtils::HomoLumoGapException& e) {
    Log::warning(e.what());
    return -1;
  }
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

void LcaoMethod::verifyPesValidity() {
  verifyChargeValidity();
  verifyMultiplicityValidity();
  verifyUnrestrictedValidity();
}

void LcaoMethod::verifyChargeValidity() {
  // Verify that the charge is not too positive (not taking away more electrons than available)
  if (molecularCharge_ > nElectronsForUnchargedSpecies_) {
    Log::error() << "The chosen molecular charge (" << molecularCharge_ << ") is too positive. Setting it to "
                 << nElectronsForUnchargedSpecies_ << ".";
    molecularCharge_ = nElectronsForUnchargedSpecies_;
  }
  else if (nElectronsForUnchargedSpecies_ - molecularCharge_ > 2 * nAOs_) {
    auto newCharge = (2 * nAOs_ - nElectronsForUnchargedSpecies_) * (-1);
    Log::error() << "Not enough orbitals to accommodate the chosen molecular charge (" << molecularCharge_
                 << "). Setting it to " << newCharge << ".";
    molecularCharge_ = newCharge;
  }
  nElectrons_ = nElectronsForUnchargedSpecies_ - molecularCharge_;
}

void LcaoMethod::verifyMultiplicityValidity() {
  if (spinMultiplicity_ > nElectrons_ + 1) {
    Log::warning() << "The chosen spin multiplicity (" << spinMultiplicity_
                   << ") is too large (not enough electrons). Setting it to " << nElectrons_ + 1 << ".";
    spinMultiplicity_ = nElectrons_ + 1;
  }
  int numberSpotsLeftForElectrons = 2 * nAOs_ - nElectrons_;
  if (spinMultiplicity_ > numberSpotsLeftForElectrons + 1) {
    Log::warning() << "The chosen spin multiplicity (" << spinMultiplicity_
                   << ") is too large (not enough orbitals). Setting it to " << numberSpotsLeftForElectrons + 1 << ".";
    spinMultiplicity_ = numberSpotsLeftForElectrons + 1;
  }
  // Check that number of electrons and spin multiplicity are compatible
  else if ((spinMultiplicity_ + nElectrons_) % 2 == 0) {
    int newMultiplicity = 1;
    if (nElectrons_ % 2 == 1)
      newMultiplicity = 2;
    Log::warning() << "The chosen spin multiplicity (" << spinMultiplicity_ << ") is not compatible with the molecular charge ("
                   << molecularCharge_ << "). Setting it to " << newMultiplicity << ".";
    spinMultiplicity_ = newMultiplicity;
  }
}

void LcaoMethod::verifyUnrestrictedValidity() {
  if (spinMultiplicity_ > 1 && !unrestrictedCalculationRunning_) {
    Log::warning() << "The chosen spin multiplicity (" << spinMultiplicity_
                   << ") requires an unrestricted calculation. Setting UHF calculation.";
    setUnrestrictedCalculation(true);
    if (!unrestrictedCalculationRunning_)
      throw UnrestrictedCalculationNotAvailableException();
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

void LcaoMethod::computeEnergyAndDerivatives(Utils::derivativeType d) {
  electronicEnergy_ = electronicPart_->calculateElectronicEnergy();
  repulsionEnergy_ = rep_->getRepulsionEnergy();
  energy_ = electronicEnergy_ + repulsionEnergy_;

  if (d == derivativeType::first)
    calculateDerivatives<Utils::derivativeType::first>(gradients_);
  else if (d == derivativeType::second_atomic)
    calculateDerivatives<Utils::derivativeType::second_atomic>(secondDerivatives_);
  else if (d == derivativeType::second_full)
    calculateDerivatives<Utils::derivativeType::second_full>(fullSecondDerivatives_);
}

template<Utils::derivativeType O>
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
  coreCharges_ = initializer_->getCoreCharges();
  unrestrictedCalculationPossible_ = initializer_->unrestrictedCalculationPossible();

  overlapCalculator_->reset();
  electronicPart_->initialize();
  rep_->initialize();

  if (!unrestrictedCalculationPossible_)
    setUnrestrictedCalculation(false);
  verifyPesValidity();
  resizeLcaoMethodMatrices();
  resizeRealTimeMethodMembers();
}

void LcaoMethod::calculateDensityIndependentQuantities(Utils::derivativeType d) {
  derivOrder order = derivOrder::zero;
  if (d == derivativeType::first)
    order = derivOrder::one;
  if (d == derivativeType::second_atomic || d == derivativeType::second_full)
    order = derivOrder::two;

  overlapCalculator_->calculateOverlap(order);
  overlapMatrix_ = overlapCalculator_->getOverlap().getMatrixXd();
  electronicPart_->calculateDensityIndependentPart(order);
  rep_->calculateRepulsion(order);
}

void LcaoMethod::assembleFockMatrix() {
  fockMatrix_ = electronicPart_->getMatrix();
}

void LcaoMethod::calculate(Utils::derivativeType d) {
  verifyPesValidity();

  calculateDensityIndependentQuantities(d);
  assembleFockMatrix();
  if (solvesOnlyOccupiedManifold()) {
    LcaoUtils::solveOccupiedRestrictedGeneralizedEigenvalueProblem(fockMatrix_, overlapMatrix_, eigenvectorMatrix_,
                                                                   singleParticleEnergies_, nElectrons_);
  }
  else {
    LcaoUtils::solveRestrictedGeneralizedEigenvalueProblem(fockMatrix_, overlapMatrix_, eigenvectorMatrix_,
                                                           singleParticleEnergies_);
  }

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
} // namespace Utils
} // namespace Scine
