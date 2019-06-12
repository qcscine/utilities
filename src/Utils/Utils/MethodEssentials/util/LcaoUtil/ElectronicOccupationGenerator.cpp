/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ElectronicOccupationGenerator.h"
#include <Utils/MethodEssentials/Methods/LCAOMethod.h>
#include <Utils/MethodEssentials/util/LcaoUtil/ElectronicOccupation.h>
#include <Utils/MethodEssentials/util/LcaoUtil/LcaoUtil.h>

namespace Scine {
namespace Utils {

namespace LcaoUtil {

void ElectronicOccupationGenerator::setMethod(LCAOMethod* method) {
  method_ = method;
}

ElectronicOccupation ElectronicOccupationGenerator::generateOccupation() {
  auto occupation = generateOccupationImpl();
  assert(occupationFulfillsRequirements(occupation));
  return occupation;
}

void ElectronicOccupationGenerator::newScfCycleStarted() {
  newScfCycleStartedImpl();
}

bool ElectronicOccupationGenerator::occupationFulfillsRequirements(const ElectronicOccupation& occupation) {
  if (method_->unrestrictedCalculationRunning()) {
    bool b1 = occupation.numberRestrictedElectrons() == 0;
    int nAlpha, nBeta;
    getNumberUnrestrictedElectrons(nAlpha, nBeta, method_->getNumberElectrons(), method_->spinMultiplicity());
    bool b2 = occupation.numberAlphaElectrons() == nAlpha;
    bool b3 = occupation.numberBetaElectrons() == nBeta;
    return b1 && b2 && b3;
  }
  else {
    bool b1 = occupation.numberAlphaElectrons() == 0;
    bool b2 = occupation.numberBetaElectrons() == 0;
    bool b3 = occupation.numberRestrictedElectrons() == method_->getNumberElectrons();
    return b1 && b2 && b3;
  }
}

} // namespace LcaoUtil
} // namespace Utils
} // namespace Scine
