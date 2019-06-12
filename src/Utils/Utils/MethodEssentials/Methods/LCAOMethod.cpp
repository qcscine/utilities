/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LCAOMethod.h"
#include "OverlapCalculator.h"
#include "StructureDependentInitializer.h"
#include <Utils/IO/Logger.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Utils/MethodEssentials/Methods/ElectronicContributionCalculator.h>
#include <Utils/MethodEssentials/Methods/RepulsionCalculator.h>
#include <Utils/MethodEssentials/util/LcaoUtil/AufbauPrincipleOccupationGenerator.h>
#include <Utils/MethodEssentials/util/LcaoUtil/DensityMatrixGenerator.h>
#include <Utils/MethodEssentials/util/LcaoUtil/ElectronicOccupationGenerator.h>
#include <Utils/MethodEssentials/util/LcaoUtil/HomoLumoGapCalculator.h>
#include <Utils/MethodEssentials/util/LcaoUtil/LcaoUtil.h>
#include <Utils/MethodEssentials/util/MatrixWithDerivatives.h>
#include <Utils/Typenames.h>
#include <iostream>

namespace Scine {
namespace Utils {

using namespace Utils::AutomaticDifferentiation;

LCAOMethod::LCAOMethod(bool unrestrictedCalculationPossible, Utils::derivOrder maximalOrder, bool orthogonalBasisSet)
  : SinglePointMethod(maximalOrder),
    nAOs_(0),
    unrestrictedCalculationPossible_(unrestrictedCalculationPossible),
    basisSetIsOrthogonal_(orthogonalBasisSet) {
  setElectronicOccupationGenerator(std::make_unique<LcaoUtil::AufbauPrincipleOccupationGenerator>());
}

LCAOMethod::~LCAOMethod() = default;

const Eigen::MatrixXd& LCAOMethod::getOverlapMatrix() const {
  return overlapMatrix_;
}

void LCAOMethod::setOverlapMatrix(const Eigen::MatrixXd& S) {
  assert(S.rows() == nAOs_ && S.cols() == nAOs_ &&
         "Dimensions of given overlap matrix do not correspond to the number of atomic orbitals.");
  overlapMatrix_ = S;
}

void LCAOMethod::setFockMatrix(SpinAdaptedMatrix F) {
  // assert(F.rows() == nAOs_ && F.cols() == nAOs_ && "Dimensions of given Fock matrix do not correspond to the number
  // of atomic orbitals."); TODO
  fockMatrix_ = std::move(F);
}

void LCAOMethod::setDensityMatrix(DensityMatrix P) {
  // assert(P.rows() == nAOs_ && P.cols() == nAOs_ && "Dimensions of given density matrix do not correspond to the
  // number of atomic orbitals."); TODO
  densityMatrix_ = std::move(P);
  if (unrestrictedCalculationRunning() && !densityMatrix_.unrestricted())
    densityMatrix_.setUnrestricted(true);
}

void LCAOMethod::setEnergyWeightedDensityMatrix(const Eigen::MatrixXd& W) {
  assert(W.rows() == nAOs_ && W.cols() == nAOs_ &&
         "Dimensions of given energy-weighted density matrix do not correspond to the number of atomic orbitals.");
  energyWeightedDensityMatrix_ = W;
}

void LCAOMethod::setMolecularOrbitals(MolecularOrbitals C) {
  // assert(C.rows() == nAOs_ && C.cols() == nAOs_ && "Dimensions of given eigenvector matrix do not correspond to the
  // number of atomic orbitals."); TODO
  eigenvectorMatrix_ = std::move(C);
}

void LCAOMethod::setUnrestrictedCalculation(bool b) {
  if (b && !unrestrictedCalculationPossible()) {
    Utils::Log::error() << "Not possible to run unrestricted calculations with this method.";
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

void LCAOMethod::setSpinMultiplicity(int s) {
  assert(s > 0 && "The spin multiplicity must be larger than zero.");
  spinMultiplicity_ = s;
}

double LCAOMethod::getHomoLumoGap() const {
  try {
    return LcaoUtil::HomoLumoGapCalculator::calculate(singleParticleEnergies_, occupation_);
  }
  catch (LcaoUtil::HomoLumoGapException& e) {
    Utils::Log::warning(e.what());
    return -1;
  }
}

void LCAOMethod::resizeLCAOMethodMatrices() {
  fockMatrix_.resize(nAOs_);
  overlapMatrix_.resize(nAOs_, nAOs_);
  densityMatrix_.resize(nAOs_);
  energyWeightedDensityMatrix_.resize(nAOs_, nAOs_);
  eigenvectorMatrix_.invalidate();
}

void LCAOMethod::setElectronicOccupationGenerator(std::unique_ptr<LcaoUtil::ElectronicOccupationGenerator>&& electronicOccupationSetter) {
  electronicOccupationGenerator_ = move(electronicOccupationSetter);
  electronicOccupationGenerator_->setMethod(this);
}

void LCAOMethod::calculateOccupation() {
  occupation_ = electronicOccupationGenerator_->generateOccupation();
}

void LCAOMethod::calculateDensity() {
  densityMatrix_ = LcaoUtil::DensityMatrixGenerator::generate(occupation_, eigenvectorMatrix_);
}

void LCAOMethod::calculateEnergyWeightedDensity() {
  energyWeightedDensityMatrix_ =
      LcaoUtil::DensityMatrixGenerator::generateEnergyWeighted(occupation_, eigenvectorMatrix_, singleParticleEnergies_);
}

void LCAOMethod::verifyPesValidity() {
  verifyChargeValidity();
  verifyMultiplicityValidity();
  verifyUnrestrictedValidity();
}

void LCAOMethod::verifyChargeValidity() {
  // Verify that the charge is not too positive (not taking away more electrons than available)
  if (molecularCharge_ > nElectronsForUnchargedSpecies_) {
    Utils::Log::error() << "The chosen molecular charge (" << molecularCharge_ << ") is too positive. Setting it to "
                        << nElectronsForUnchargedSpecies_ << ".";
    molecularCharge_ = nElectronsForUnchargedSpecies_;
  }
  else if (nElectronsForUnchargedSpecies_ - molecularCharge_ > 2 * nAOs_) {
    auto newCharge = (2 * nAOs_ - nElectronsForUnchargedSpecies_) * (-1);
    Utils::Log::error() << "Not enough orbitals to accommodate the chosen molecular charge (" << molecularCharge_
                        << "). Setting it to " << newCharge << ".";
    molecularCharge_ = newCharge;
  }
  nElectrons_ = nElectronsForUnchargedSpecies_ - molecularCharge_;
}

void LCAOMethod::verifyMultiplicityValidity() {
  if (spinMultiplicity_ > nElectrons_ + 1) {
    Utils::Log::warning() << "The chosen spin multiplicity (" << spinMultiplicity_
                          << ") is too large (not enough electrons). Setting it to " << nElectrons_ + 1 << ".";
    spinMultiplicity_ = nElectrons_ + 1;
  }
  int numberSpotsLeftForElectrons = 2 * nAOs_ - nElectrons_;
  if (spinMultiplicity_ > numberSpotsLeftForElectrons + 1) {
    Utils::Log::warning() << "The chosen spin multiplicity (" << spinMultiplicity_
                          << ") is too large (not enough orbitals). Setting it to " << numberSpotsLeftForElectrons + 1
                          << ".";
    spinMultiplicity_ = numberSpotsLeftForElectrons + 1;
  }
  // Check that number of electrons and spin multiplicity are compatible
  else if ((spinMultiplicity_ + nElectrons_) % 2 == 0) {
    int newMultiplicity = 1;
    if (nElectrons_ % 2 == 1)
      newMultiplicity = 2;
    Utils::Log::warning() << "The chosen spin multiplicity (" << spinMultiplicity_
                          << ") is not compatible with the molecular charge (" << molecularCharge_
                          << "). Setting it to " << newMultiplicity << ".";
    spinMultiplicity_ = newMultiplicity;
  }
}

void LCAOMethod::verifyUnrestrictedValidity() {
  if (spinMultiplicity_ > 1 && !unrestrictedCalculationRunning_) {
    Utils::Log::warning() << "The chosen spin multiplicity (" << spinMultiplicity_
                          << ") requires an unrestricted calculation. Setting UHF calculation.";
    setUnrestrictedCalculation(true);
  }
}

void LCAOMethod::calculateBondOrderMatrix() {
  if (basisSetIsOrthogonal())
    LcaoUtil::calculateOrthonormalBondOrderMatrix(bondOrders_, densityMatrix_, aoIndexes_);
  else
    LcaoUtil::calculateBondOrderMatrix(bondOrders_, densityMatrix_, overlapMatrix_, aoIndexes_);
}

void LCAOMethod::calculateAtomicCharges() {
  if (basisSetIsOrthogonal())
    LcaoUtil::calculateOrthonormalAtomicCharges(atomicCharges_, coreCharges_, densityMatrix_, aoIndexes_);
  else
    LcaoUtil::calculateMullikenAtomicCharges(atomicCharges_, coreCharges_, densityMatrix_, overlapMatrix_, aoIndexes_);
}

void LCAOMethod::computeEnergyAndDerivatives(Utils::derivativeType d) {
  electronicEnergy_ = electronicPart_->calculateElectronicEnergy();
  repulsionEnergy_ = rep_->getRepulsionEnergy();
  energy_ = electronicEnergy_ + repulsionEnergy_;

  if (d == Utils::derivativeType::first)
    calculateDerivatives<Utils::derivativeType::first>(gradients_);
  else if (d == Utils::derivativeType::second_atomic)
    calculateDerivatives<Utils::derivativeType::second_atomic>(secondDerivatives_);
  else if (d == Utils::derivativeType::second_full)
    calculateDerivatives<Utils::derivativeType::second_full>(fullSecondDerivatives_);
}

template<Utils::derivativeType O>
void LCAOMethod::calculateDerivatives(DerivativeContainerType<O>& derivatives) {
  derivatives.setZero();

  rep_->addRepulsionDerivatives(derivatives);
  electronicPart_->addDerivatives(derivatives);
}

void LCAOMethod::initialize() {
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
  resizeLCAOMethodMatrices();
  resizeRealTimeMethodMembers();
}

void LCAOMethod::calculateDensityIndependentQuantities(Utils::derivativeType d) {
  Utils::derivOrder order = Utils::derivOrder::zero;
  if (d == Utils::derivativeType::first)
    order = Utils::derivOrder::one;
  if (d == Utils::derivativeType::second_atomic || d == Utils::derivativeType::second_full)
    order = Utils::derivOrder::two;

  overlapCalculator_->calculateOverlap(order);
  overlapMatrix_ = overlapCalculator_->getOverlap().getMatrixXd();
  electronicPart_->calculateDensityIndependentPart(order);
  rep_->calculateRepulsion(order);
}

void LCAOMethod::assembleFockMatrix() {
  fockMatrix_ = electronicPart_->getMatrix();
}

// TODO: Davidson for large systems
void LCAOMethod::calculate(Utils::derivativeType d) {
  verifyPesValidity();

  calculateDensityIndependentQuantities(d);
  assembleFockMatrix();
  LcaoUtil::solveRestrictedGeneralizedEigenvalueProblem(fockMatrix_, overlapMatrix_, eigenvectorMatrix_, singleParticleEnergies_);

  calculateOccupationAndDensity();
  calculateBondOrderMatrix();
  calculateAtomicCharges();

  computeEnergyAndDerivatives(d);
}

void LCAOMethod::calculateOccupationAndDensity() {
  calculateOccupation();
  calculateDensity();
  if (!basisSetIsOrthogonal()) {
    calculateEnergyWeightedDensity();
  }
}
} // namespace Utils
} // namespace Scine
