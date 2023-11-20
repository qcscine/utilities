/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AufbauPrincipleOccupationGenerator.h"
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>
#include <Utils/Scf/MethodInterfaces/ScfMethod.h>

namespace Scine {
namespace Utils {

namespace LcaoUtils {

ElectronicOccupation AufbauPrincipleOccupationGenerator::generateOccupationImpl() {
  ElectronicOccupation occupation;

  if (method_->unrestrictedCalculationRunning()) {
    int nAlpha, nBeta;
    getNumberUnrestrictedElectrons(nAlpha, nBeta, method_->getNumberElectrons(), method_->spinMultiplicity());
    occupation.fillLowestUnrestrictedOrbitals(nAlpha, nBeta);
  }
  else {
    occupation.fillLowestRestrictedOrbitalsWithElectrons(method_->getNumberElectrons());
  }

  return occupation;
}

} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine
