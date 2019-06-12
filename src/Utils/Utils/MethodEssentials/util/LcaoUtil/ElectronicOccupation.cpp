/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ElectronicOccupation.h"
#include <cassert>
#include <numeric>

namespace Scine {
namespace Utils {

namespace LcaoUtil {

int ElectronicOccupation::numberRestrictedElectrons() const {
  return nRestrictedElectrons_;
}

int ElectronicOccupation::numberOccupiedRestrictedOrbitals() const {
  return (nRestrictedElectrons_ + 1) / 2;
}

int ElectronicOccupation::numberAlphaElectrons() const {
  return nAlphaElectrons_;
}

int ElectronicOccupation::numberBetaElectrons() const {
  return nBetaElectrons_;
}

const std::vector<int>& ElectronicOccupation::getFilledRestrictedOrbitals() const {
  if (!specifiedWithOrbitals_)
    filledRestrictedOrbitals_ = generateOccupiedArrayFromNumberOfOccupiedOrbitals(nRestrictedElectrons_ / 2);
  return filledRestrictedOrbitals_;
}

const std::vector<int>& ElectronicOccupation::getFilledAlphaOrbitals() const {
  if (!specifiedWithOrbitals_)
    filledAlphaOrbitals_ = generateOccupiedArrayFromNumberOfOccupiedOrbitals(nAlphaElectrons_);
  return filledAlphaOrbitals_;
}

const std::vector<int>& ElectronicOccupation::getFilledBetaOrbitals() const {
  if (!specifiedWithOrbitals_)
    filledBetaOrbitals_ = generateOccupiedArrayFromNumberOfOccupiedOrbitals(nBetaElectrons_);
  return filledBetaOrbitals_;
}

bool ElectronicOccupation::isFilledUpFromTheBottom() const {
  if (!specifiedWithNumberElectrons_)
    checkWhetherFilledFromBottom();
  return specifiedWithNumberElectrons_;
}

bool ElectronicOccupation::hasUnpairedRHFElectron() const {
  return hasUnpairedElectron_;
}

bool ElectronicOccupation::isUnrestricted() const {
  return !isRestricted();
}

bool ElectronicOccupation::isRestricted() const {
  return restricted_;
}

void ElectronicOccupation::makeUnrestricted() {
  if (isUnrestricted())
    return; // Nothing to be done

  auto restrictedOrbitals = getFilledRestrictedOrbitals();
  fillSpecifiedUnrestrictedOrbitals(restrictedOrbitals, restrictedOrbitals);
}

ElectronicOccupation ElectronicOccupation::toUnrestricted() const {
  ElectronicOccupation eo = *this;
  eo.makeUnrestricted();
  return eo;
}

void ElectronicOccupation::fillLowestRestrictedOrbitalsWithElectrons(int nElectrons) {
  assert(nElectrons >= 0);

  reset();

  restricted_ = true;
  specifiedWithNumberElectrons_ = true;
  nRestrictedElectrons_ = nElectrons;
  if (nRestrictedElectrons_ % 2 == 1)
    hasUnpairedElectron_ = true;
}

void ElectronicOccupation::fillLowestUnrestrictedOrbitals(int nAlphaElectrons, int nBetaElectrons) {
  assert(nAlphaElectrons >= 0 && nBetaElectrons >= 0);

  reset();

  restricted_ = false;
  specifiedWithNumberElectrons_ = true;
  nAlphaElectrons_ = nAlphaElectrons;
  nBetaElectrons_ = nBetaElectrons;
}

void ElectronicOccupation::fillSpecifiedRestrictedOrbitals(std::vector<int> occupiedRestrictedOrbitals) {
  reset();
  restricted_ = true;
  specifiedWithOrbitals_ = true;
  filledRestrictedOrbitals_ = std::move(occupiedRestrictedOrbitals);
  nRestrictedElectrons_ = 2 * static_cast<int>(filledRestrictedOrbitals_.size());
}

void ElectronicOccupation::fillSpecifiedUnrestrictedOrbitals(std::vector<int> occupiedAlphaOrbitals,
                                                             std::vector<int> occupiedBetaOrbitals) {
  reset();
  restricted_ = false;
  specifiedWithOrbitals_ = true;
  filledAlphaOrbitals_ = std::move(occupiedAlphaOrbitals);
  filledBetaOrbitals_ = std::move(occupiedBetaOrbitals);
  nAlphaElectrons_ = static_cast<int>(filledAlphaOrbitals_.size());
  nBetaElectrons_ = static_cast<int>(filledBetaOrbitals_.size());
}

void ElectronicOccupation::reset() {
  *this = ElectronicOccupation{};
}

std::vector<int> ElectronicOccupation::generateOccupiedArrayFromNumberOfOccupiedOrbitals(int nOrbitals) const {
  std::vector<int> orbitals(static_cast<std::size_t>(nOrbitals));
  // std::iota generates a range of consecutive numbers
  std::iota(orbitals.begin(), orbitals.end(), 0);
  specifiedWithNumberElectrons_ = true;
  return orbitals;
}

void ElectronicOccupation::checkWhetherFilledFromBottom() const {
  specifiedWithNumberElectrons_ = true;
  for (int i = 0; i < static_cast<int>(filledRestrictedOrbitals_.size()); ++i)
    if (filledRestrictedOrbitals_[i] != i)
      specifiedWithNumberElectrons_ = false;
  for (int i = 0; i < static_cast<int>(filledAlphaOrbitals_.size()); ++i)
    if (filledAlphaOrbitals_[i] != i)
      specifiedWithNumberElectrons_ = false;
  for (int i = 0; i < static_cast<int>(filledBetaOrbitals_.size()); ++i)
    if (filledBetaOrbitals_[i] != i)
      specifiedWithNumberElectrons_ = false;
}

} // namespace LcaoUtil
} // namespace Utils
} // namespace Scine
