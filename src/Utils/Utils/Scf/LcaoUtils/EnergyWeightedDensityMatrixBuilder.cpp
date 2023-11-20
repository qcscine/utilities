/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "EnergyWeightedDensityMatrixBuilder.h"
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>
#include <cassert>
#include <cmath>

namespace Scine {
namespace Utils {

namespace LcaoUtils {

EnergyWeightedDensityMatrixBuilder::EnergyWeightedDensityMatrixBuilder(const MolecularOrbitals& coefficientMatrix,
                                                                       const SingleParticleEnergies& orbitalEnergies)
  : coefficientMatrix_(coefficientMatrix), orbitalEnergies_(orbitalEnergies) {
  assert(coefficientMatrix.isValid());
}

Eigen::MatrixXd EnergyWeightedDensityMatrixBuilder::generateRestrictedForNumberElectrons(int nElectrons) const {
  assert(nElectrons >= 0);
  assert(coefficientMatrix_.isRestricted());
  const auto& C = coefficientMatrix_.restrictedMatrix();
  auto nAOs = C.rows();
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(nAOs, nAOs);

  for (int i = 0; i < nElectrons / 2; ++i) {
    m += 2 * calculateSingleOrbitalEWDensity(C.col(i), orbitalEnergies_.getRestrictedLevelEnergy(i));
  }

  // if odd number of electrons
  if ((nElectrons % 2) != 0) {
    m += calculateSingleOrbitalEWDensity(C.col(nElectrons / 2), orbitalEnergies_.getRestrictedLevelEnergy(nElectrons / 2));
  }

  return m;
}

Eigen::MatrixXd
EnergyWeightedDensityMatrixBuilder::generateUnrestrictedForNumberElectronsAndMultiplicity(int nElectrons,
                                                                                          int spinMultiplicity) const {
  int nAlpha, nBeta;
  LcaoUtils::getNumberUnrestrictedElectrons(nAlpha, nBeta, nElectrons, spinMultiplicity);
  return generateUnrestrictedForNumberAlphaAndBetaElectrons(nAlpha, nBeta);
}

Eigen::MatrixXd EnergyWeightedDensityMatrixBuilder::generateUnrestrictedForNumberAlphaAndBetaElectrons(int nAlpha,
                                                                                                       int nBeta) const {
  assert(coefficientMatrix_.isUnrestricted());
  const auto& cA = coefficientMatrix_.alphaMatrix();
  const auto& cB = coefficientMatrix_.betaMatrix();
  auto nAOs = cA.rows();
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(nAOs, nAOs);

  for (int i = 0; i < nAlpha; ++i) {
    m += calculateSingleOrbitalEWDensity(cA.col(i), orbitalEnergies_.getAlphaLevelEnergy(i));
  }
  for (int i = 0; i < nBeta; ++i) {
    m += calculateSingleOrbitalEWDensity(cB.col(i), orbitalEnergies_.getBetaLevelEnergy(i));
  }

  return m;
}

Eigen::MatrixXd
EnergyWeightedDensityMatrixBuilder::generateRestrictedForSpecifiedOrbitals(const std::vector<int>& doublyOccupiedOrbitals) const {
  assert(coefficientMatrix_.isRestricted());
  const auto& C = coefficientMatrix_.restrictedMatrix();
  auto nAOs = C.rows();
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(nAOs, nAOs);

  for (auto o : doublyOccupiedOrbitals) {
    assert(o < nAOs && "Orbital index larger than the number of atomic orbitals.");
    m += 2 * calculateSingleOrbitalEWDensity(C.col(o), orbitalEnergies_.getRestrictedLevelEnergy(o));
  }

  return m;
}

Eigen::MatrixXd
EnergyWeightedDensityMatrixBuilder::generateUnrestrictedForSpecifiedOrbitals(const std::vector<int>& alphaOrbitals,
                                                                             const std::vector<int>& betaOrbitals) const {
  assert(coefficientMatrix_.isUnrestricted());
  const auto& Ca = coefficientMatrix_.alphaMatrix();
  const auto& Cb = coefficientMatrix_.betaMatrix();
  auto nAOs = Ca.rows();
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(nAOs, nAOs);

  for (auto o : alphaOrbitals) {
    assert(o < nAOs && "Orbital index larger than the number of atomic orbitals.");
    m += calculateSingleOrbitalEWDensity(Ca.col(o), orbitalEnergies_.getAlphaLevelEnergy(o));
  }
  for (auto o : betaOrbitals) {
    assert(o < nAOs && "Orbital index larger than the number of atomic orbitals.");
    m += calculateSingleOrbitalEWDensity(Cb.col(o), orbitalEnergies_.getBetaLevelEnergy(o));
  }

  return m;
}

Eigen::MatrixXd EnergyWeightedDensityMatrixBuilder::calculateSingleOrbitalEWDensity(const Eigen::VectorXd& eigenvector,
                                                                                    double orbitalEnergy) {
  return orbitalEnergy * eigenvector * eigenvector.transpose();
}

} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine
