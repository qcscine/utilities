/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ElectronicOccupationGenerator.h"
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>
#include <Utils/Scf/MethodInterfaces/LcaoMethod.h>

namespace Scine {
namespace Utils {

namespace LcaoUtils {

void ElectronicOccupationGenerator::setMethod(LcaoMethod* method) {
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

  bool b1 = occupation.numberAlphaElectrons() == 0;
  bool b2 = occupation.numberBetaElectrons() == 0;
  bool b3 = occupation.numberRestrictedElectrons() == method_->getNumberElectrons();
  return b1 && b2 && b3;
}

} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine
