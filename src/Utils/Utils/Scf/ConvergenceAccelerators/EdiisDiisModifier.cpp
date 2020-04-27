/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "EdiisDiisModifier.h"
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>
#include <Utils/Scf/MethodInterfaces/ScfMethod.h>

namespace Scine {
namespace Utils {

EdiisDiisModifier::EdiisDiisModifier() {
  setSpaceSize(6);
}

void EdiisDiisModifier::setSpaceSize(unsigned n) {
  diis_.setSubspaceSize(n);
  ediis_.setSubspaceSize(n);
}

void EdiisDiisModifier::onOverlapCalculated() {
  if (!initialized) {
    initialize();
    initialized = true;
  }

  diis_.setNAOs(m->getNumberAtomicOrbitals());
  ediis_.setNAOs(m->getNumberAtomicOrbitals());
  ediis_.restart();
  diis_.setOverlapMatrix(m->getOverlapMatrix());

  if (m->unrestrictedCalculationRunning()) {
    ediis_.setUnrestricted(true);
    diis_.setUnrestricted(true);
  }
  else {
    diis_.setUnrestricted(false);
    ediis_.setUnrestricted(false);
  }
}

void EdiisDiisModifier::onFockCalculated() {
  if (m->unrestrictedCalculationRunning()) {
    if (!sameNumberOfElectronsInMethodAndInDensityMatrix())
      return;
  }
  m->computeEnergyAndDerivatives(Utils::derivativeType::none);
  ediis_.addMatrices(m->getEnergy(), m->getFockMatrix(), m->getDensityMatrix());
  diis_.addMatrices(m->getFockMatrix(), m->getDensityMatrix());

  m->setFockMatrix(getCombinedFockMatrix());
}

void EdiisDiisModifier::initialize() {
  if (m->basisSetIsOrthogonal())
    setOrthogonal(true);
}

void EdiisDiisModifier::setOrthogonal(bool o) {
  diis_.setOrthogonal(o);
}

bool EdiisDiisModifier::sameNumberOfElectronsInMethodAndInDensityMatrix() {
  int nAlphaInMethod, nBetaInMethod;
  LcaoUtils::getNumberUnrestrictedElectrons(nAlphaInMethod, nBetaInMethod, m->getNumberElectrons(), m->spinMultiplicity());

  int nAlphaInMatrix = m->getDensityMatrix().numberElectronsInAlphaMatrix();
  int nBetaInMatrix = m->getDensityMatrix().numberElectronsInBetaMatrix();

  return nAlphaInMethod == nAlphaInMatrix && nBetaInMethod == nBetaInMatrix;
}

SpinAdaptedMatrix EdiisDiisModifier::getCombinedFockMatrix() {
  double errMax = diis_.getMaxError();
  double errMin = diis_.getMinError();
  double errLast = diis_.getLastError();

  if (errMax > 1e-1 || errLast > 1.1 * errMin)
    return ediis_.getMixedFockMatrix();
  if (errMax < 1e-4)
    return diis_.getMixedFockMatrix();

  // Else: linear combination
  return mixedFockMatrix(errMax);
}

SpinAdaptedMatrix EdiisDiisModifier::mixedFockMatrix(double errMax) {
  double cEdiis = 10 * errMax;
  double cDiis = (1. - 10 * errMax);
  auto fEdiis = ediis_.getMixedFockMatrix();
  auto fDiis = diis_.getMixedFockMatrix();
  if (m->unrestrictedCalculationRunning())
    return SpinAdaptedMatrix::createUnrestricted(fDiis.alphaMatrix() * cDiis + fEdiis.alphaMatrix() * cEdiis,
                                                 fDiis.betaMatrix() * cDiis + fEdiis.betaMatrix() * cEdiis);
  else
    return SpinAdaptedMatrix::createRestricted(fDiis.restrictedMatrix() * cDiis + fEdiis.restrictedMatrix() * cEdiis);
}
} // namespace Utils
} // namespace Scine
