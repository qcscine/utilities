/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "EdiisModifier.h"
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>
#include <Utils/Scf/MethodInterfaces/ScfMethod.h>
#include <iostream>

namespace Scine {
namespace Utils {

void EdiisModifier::onOverlapCalculated() {
  mixer_.setNAOs(m->getNumberAtomicOrbitals());
  mixer_.restart();

  if (m->unrestrictedCalculationRunning()) {
    mixer_.setUnrestricted(true);
  }
}

void EdiisModifier::onFockCalculated() {
  if (!sameNumberOfElectronsInMethodAndInDensityMatrix()) {
    return;
  }
  m->computeEnergyAndDerivatives(Utils::Derivative::None);
  mixer_.addMatrices(m->getEnergy(), m->getFockMatrix(), m->getDensityMatrix());
  m->setFockMatrix(mixer_.getMixedFockMatrix());
}

void EdiisModifier::setSpaceSize(int n) {
  mixer_.setSubspaceSize(n);
}

bool EdiisModifier::sameNumberOfElectronsInMethodAndInDensityMatrix() {
  int nAlphaInMethod, nBetaInMethod;
  LcaoUtils::getNumberUnrestrictedElectrons(nAlphaInMethod, nBetaInMethod, m->getNumberElectrons(), m->spinMultiplicity());

  int nAlphaInMatrix = m->getDensityMatrix().numberElectronsInAlphaMatrix();
  int nBetaInMatrix = m->getDensityMatrix().numberElectronsInBetaMatrix();

  return nAlphaInMethod == nAlphaInMatrix && nBetaInMethod == nBetaInMatrix;
}
} // namespace Utils
} // namespace Scine
