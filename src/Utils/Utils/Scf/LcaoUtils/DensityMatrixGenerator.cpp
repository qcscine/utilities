/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DensityMatrixGenerator.h"
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/Scf/LcaoUtils/DensityMatrixBuilder.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Scf/LcaoUtils/EnergyWeightedDensityMatrixBuilder.h>
#include <cassert>

namespace Scine {
namespace Utils {

namespace LcaoUtils {

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

} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine
