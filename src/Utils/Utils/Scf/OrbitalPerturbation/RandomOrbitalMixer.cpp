/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "RandomOrbitalMixer.h"
#include "UniqueRandomNumbersGenerator.h"
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <random>

namespace Scine {
namespace Utils {

namespace OrbitalPerturbation {

std::mt19937& RandomOrbitalMixer::getRandomNumberGenerator() {
  static std::mt19937 rdGen; // NB: if non-static, the same orbital mixes will always be generated.
  return rdGen;
}

RandomOrbitalMixer::RandomOrbitalMixer(MolecularOrbitals& orbitals, int nRestrictedElectrons)
  : orbitals_(orbitals),
    alphaHomo_(nRestrictedElectrons / 2 - 1),
    betaHomo_(nRestrictedElectrons / 2 - 1),
    nOrbitals_(orbitals.numberOrbitals()) {
  assert(orbitals.isRestricted());
}

RandomOrbitalMixer::RandomOrbitalMixer(MolecularOrbitals& orbitals, int nAlphaElectrons, int nBetaElectrons)
  : orbitals_(orbitals), alphaHomo_(nAlphaElectrons - 1), betaHomo_(nBetaElectrons - 1), nOrbitals_(orbitals.numberOrbitals()) {
  if (!orbitals.isUnrestricted()) {
    orbitals.makeUnrestricted();
  }
}

void RandomOrbitalMixer::mix(Core::Log& log) {
  checkValidNumberOfMixes();

  if (invalidMolecularOrbitals(orbitals_)) {
    throw std::runtime_error("Cannot mix the molecular orbitals since they are invalid.");
  }

  if (orbitals_.isUnrestricted()) {
    auto alphaMixes = calculateMixes(alphaHomo_);
    auto betaMixes = calculateMixes(betaHomo_);
    LcaoUtils::MolecularOrbitalsManipulation::mixUnrestricted(orbitals_, alphaMixes, betaMixes);
  }
  else {
    log.warning
        << "WARNING: Restricted orbitals detected. Consider an unrestricted calculation for effective orbital mixing."
        << Core::Log::endl;
    auto restrictedMixes = calculateMixes(alphaHomo_);
    LcaoUtils::MolecularOrbitalsManipulation::mixRestricted(orbitals_, restrictedMixes);
  }
}

bool RandomOrbitalMixer::invalidMolecularOrbitals(const MolecularOrbitals& mo) {
  return !mo.isValid();
}

std::vector<LcaoUtils::MolecularOrbitalsManipulation::Mix> RandomOrbitalMixer::calculateMixes(int homoIndex) const {
  auto selectedOccupied = getRandomOccupiedOrbitals(homoIndex);
  auto selectedVirtual = getRandomVirtualOrbitals(homoIndex);
  auto mixes = createMixes(selectedOccupied, selectedVirtual);
  return mixes;
}

std::vector<int> RandomOrbitalMixer::getRandomOccupiedOrbitals(int homoIndex) const {
  int minOccupied = calculateMinimalOccupiedIndex(homoIndex);
  // i.e. if there are no electrons at all
  if (homoIndex < 0) {
    return {};
  }
  return selectUniqueRandomNumbers(minOccupied, homoIndex);
}

std::vector<int> RandomOrbitalMixer::getRandomVirtualOrbitals(int homoIndex) const {
  int maxVirtual = calculateMaximalVirtualIndex(homoIndex);
  if (homoIndex >= maxVirtual) { // Avoids problems when there is no virtual orbital
    return {};
  }
  return selectUniqueRandomNumbers(homoIndex + 1, maxVirtual);
}

std::vector<int> RandomOrbitalMixer::selectUniqueRandomNumbers(int min, int max) const {
  assert(min <= max);
  UniqueRandomNumbersGenerator<int> randomGenerator;
  randomGenerator.setRange(min, max);
  return randomGenerator.generate(getRandomNumberGenerator(), static_cast<unsigned>(numberMixes_));
}

std::vector<LcaoUtils::MolecularOrbitalsManipulation::Mix>
RandomOrbitalMixer::createMixes(const std::vector<int>& occ, const std::vector<int>& virt) const {
  std::vector<LcaoUtils::MolecularOrbitalsManipulation::Mix> mixes;

  std::uniform_real_distribution<double> distr(minAngle_, maxAngle_);
  for (int i = 0; i < numberMixes_; ++i) {
    double randomNumber = distr(getRandomNumberGenerator());
    mixes.emplace_back(occ[i], virt[i], randomNumber);
  }

  return mixes;
}

int RandomOrbitalMixer::calculateMinimalOccupiedIndex(int homoIndex) const {
  if (considerAllOrbitals_) {
    return 0;
  }

  if (numberOrbitalsToConsider_ > homoIndex + 1) {
    return 0;
  }

  return homoIndex + 1 - numberOrbitalsToConsider_;
}

int RandomOrbitalMixer::calculateMaximalVirtualIndex(int homoIndex) const {
  if (considerAllOrbitals_) {
    return nOrbitals_ - 1;
  }

  if (homoIndex + numberOrbitalsToConsider_ > nOrbitals_ - 1) {
    return nOrbitals_ - 1;
  }

  return homoIndex + numberOrbitalsToConsider_;
}

void RandomOrbitalMixer::checkValidNumberOfMixes() {
  int nOccupied = std::min(alphaHomo_, betaHomo_) + 1;
  int nVirtual = nOrbitals_ - (std::max(alphaHomo_, betaHomo_) + 1);

  int nPossibleMixes = std::min(nOccupied, nVirtual);
  if (numberMixes_ > nPossibleMixes) {
    numberMixes_ = nPossibleMixes;
  }
}

void RandomOrbitalMixer::considerOnlyOrbitalsCloseToFrontierOrbitals(int numberToConsider) {
  considerAllOrbitals_ = false;
  numberOrbitalsToConsider_ = numberToConsider;
}

void RandomOrbitalMixer::considerAllOrbitals() {
  considerAllOrbitals_ = true;
}

} // namespace OrbitalPerturbation
} // namespace Utils
} // namespace Scine
