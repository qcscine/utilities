/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "HomoLumoGapCalculator.h"
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <cassert>

namespace Scine {
namespace Utils {

namespace LcaoUtils {

double HomoLumoGapCalculator::calculate(const SingleParticleEnergies& energies, const ElectronicOccupation& occupation) {
  assert(energies.isRestricted() == occupation.isRestricted());
  if (energies.isRestricted())
    return calculateRestricted(energies, occupation);
  else
    return calculateUnrestricted(energies, occupation);
}

double HomoLumoGapCalculator::calculateRestricted(const SingleParticleEnergies& energies, const ElectronicOccupation& occupation) {
  auto nRestrictedOrbitals = energies.getRestrictedNLevels();
  auto nElectrons = occupation.numberRestrictedElectrons();
  int homoIndex = (nElectrons - 1) / 2;
  int lumoIndex = homoIndex + 1;

  if (nElectrons == 0)
    throw HomoLumoGapException("The HOMO-LUMO gap cannot be calculated because there are no electrons in the system.");
  if (lumoIndex >= nRestrictedOrbitals)
    throw HomoLumoGapException("The HOMO-LUMO gap cannot be calculated because there are no empty orbitals.");

  return energies.getRestrictedLevelEnergy(lumoIndex) - energies.getRestrictedLevelEnergy(homoIndex);
}

double HomoLumoGapCalculator::calculateUnrestricted(const SingleParticleEnergies& energies,
                                                    const ElectronicOccupation& occupation) {
  auto nOrbitals = energies.getUnrestrictedNLevels();
  auto nAlphaElectrons = occupation.numberAlphaElectrons();
  auto nBetaElectrons = occupation.numberBetaElectrons();

  if (nAlphaElectrons + nBetaElectrons == 0)
    throw HomoLumoGapException("The HOMO-LUMO gap cannot be calculated because there are no electrons in the system.");

  int homoAlphaIndex = nAlphaElectrons - 1;
  int lumoAlphaIndex = nAlphaElectrons;
  int homoBetaIndex = nBetaElectrons - 1;
  int lumoBetaIndex = nBetaElectrons;
  if (lumoAlphaIndex >= nOrbitals && lumoBetaIndex >= nOrbitals)
    throw HomoLumoGapException("The HOMO-LUMO gap cannot be calculated because there are no empty orbitals.");

  double homoAlpha = std::numeric_limits<double>::min();
  double homoBeta = std::numeric_limits<double>::min();
  double lumoAlpha = std::numeric_limits<double>::max();
  double lumoBeta = std::numeric_limits<double>::max();
  if (lumoAlphaIndex < nOrbitals)
    lumoAlpha = energies.getAlphaLevelEnergy(lumoAlphaIndex);
  if (lumoBetaIndex < nOrbitals)
    lumoBeta = energies.getBetaLevelEnergy(lumoBetaIndex);
  if (homoAlphaIndex >= 0)
    homoAlpha = energies.getAlphaLevelEnergy(homoAlphaIndex);
  if (homoBetaIndex >= 0)
    homoBeta = energies.getBetaLevelEnergy(homoBetaIndex);

  double homo = std::max(homoAlpha, homoBeta);
  double lumo = std::min(lumoAlpha, lumoBeta);

  return lumo - homo;
}

} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine
