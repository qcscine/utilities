/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DensityMatrixGenerator.h"
#include <Utils/MethodEssentials/util/DensityMatrix.h>
#include <Utils/MethodEssentials/util/LcaoUtil/DensityMatrixBuilder.h>
#include <Utils/MethodEssentials/util/LcaoUtil/ElectronicOccupation.h>
#include <Utils/MethodEssentials/util/LcaoUtil/EnergyWeightedDensityMatrixBuilder.h>
#include <Utils/MethodEssentials/util/MolecularOrbitals.h>
#include <cassert>

namespace Scine {
namespace Utils {

namespace LcaoUtil {

DensityMatrix DensityMatrixGenerator::generate(const ElectronicOccupation& occupation,
                                               const MolecularOrbitals& molecularOrbitals) {
  assert(molecularOrbitals.isValid());
  DensityMatrixBuilder builder(molecularOrbitals);
  if (occupation.isUnrestricted()) {
    if (occupation.isFilledUpFromTheBottom())
      return builder.generateUnrestrictedForNumberAlphaAndBetaElectrons(occupation.numberAlphaElectrons(),
                                                                        occupation.numberBetaElectrons());
    else
      return builder.generateUnrestrictedForSpecifiedOrbitals(occupation.getFilledAlphaOrbitals(),
                                                              occupation.getFilledBetaOrbitals());
  }
  else {
    if (occupation.isFilledUpFromTheBottom())
      return builder.generateRestrictedForNumberElectrons(occupation.numberRestrictedElectrons());
    else
      return builder.generateRestrictedForSpecifiedOrbitals(occupation.getFilledRestrictedOrbitals());
  }
}

Eigen::MatrixXd DensityMatrixGenerator::generateEnergyWeighted(const ElectronicOccupation& occupation,
                                                               const MolecularOrbitals& molecularOrbitals,
                                                               const SingleParticleEnergies& orbitalEnergies) {
  assert(molecularOrbitals.isValid());
  EnergyWeightedDensityMatrixBuilder builder(molecularOrbitals, orbitalEnergies);
  if (occupation.isUnrestricted()) {
    if (occupation.isFilledUpFromTheBottom())
      return builder.generateUnrestrictedForNumberAlphaAndBetaElectrons(occupation.numberAlphaElectrons(),
                                                                        occupation.numberBetaElectrons());
    else
      return builder.generateUnrestrictedForSpecifiedOrbitals(occupation.getFilledAlphaOrbitals(),
                                                              occupation.getFilledBetaOrbitals());
  }
  else {
    if (occupation.isFilledUpFromTheBottom())
      return builder.generateRestrictedForNumberElectrons(occupation.numberRestrictedElectrons());
    else
      return builder.generateRestrictedForSpecifiedOrbitals(occupation.getFilledRestrictedOrbitals());
  }
}

} // namespace LcaoUtil
} // namespace Utils
} // namespace Scine
