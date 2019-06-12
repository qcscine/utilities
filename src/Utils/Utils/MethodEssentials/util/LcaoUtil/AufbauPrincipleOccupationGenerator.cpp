/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AufbauPrincipleOccupationGenerator.h"
#include <Utils/MethodEssentials/Methods/SCFMethod.h>
#include <Utils/MethodEssentials/util/LcaoUtil/ElectronicOccupation.h>
#include <Utils/MethodEssentials/util/LcaoUtil/LcaoUtil.h>

namespace Scine {
namespace Utils {

namespace LcaoUtil {

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

} // namespace LcaoUtil
} // namespace Utils
} // namespace Scine
