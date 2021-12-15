/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DensityMatrixBuilder.h"
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>
#include <Eigen/Core>
#include <cassert>
#include <cmath>

namespace Scine {
namespace Utils {

namespace LcaoUtils {

DensityMatrixBuilder::DensityMatrixBuilder(const MolecularOrbitals& coefficientMatrix)
  : coefficientMatrix_(coefficientMatrix) {
  assert(coefficientMatrix.isValid());
}

DensityMatrix DensityMatrixBuilder::generateRestrictedForNumberElectrons(int nElectrons) const {
  assert(nElectrons >= 0);
  assert(coefficientMatrix_.isRestricted());
  const auto& C = coefficientMatrix_.restrictedMatrix();
  Eigen::MatrixXd P = 2 * calculateDensityMatrix(C, nElectrons / 2);

  // if odd number of electrons
  if ((nElectrons % 2) != 0) {
    P += calculateSingleOrbitalDensity(C.col(nElectrons / 2));
  }

  DensityMatrix densityMatrix;
  densityMatrix.setDensity(std::move(P), nElectrons);
  return densityMatrix;
}

Eigen::MatrixXd DensityMatrixBuilder::calculateDensityMatrix(const Eigen::MatrixXd& coefficientMatrix, int nOccupiedLevels) {
  auto nAOs = coefficientMatrix.cols();
  assert(nAOs >= nOccupiedLevels && "More electrons than atomic orbitals.");

  Eigen::MatrixXd P = calculateBlockOrbitalDensity(coefficientMatrix.block(0, 0, nAOs, nOccupiedLevels));
  return P;
}

DensityMatrix DensityMatrixBuilder::generateUnrestrictedForNumberElectronsAndMultiplicity(int nElectrons,
                                                                                          int spinMultiplicity) const {
  assert(coefficientMatrix_.isUnrestricted());
  int nAlpha, nBeta;
  LcaoUtils::getNumberUnrestrictedElectrons(nAlpha, nBeta, nElectrons, spinMultiplicity);
  return generateUnrestrictedForNumberAlphaAndBetaElectrons(nAlpha, nBeta);
}

DensityMatrix DensityMatrixBuilder::generateUnrestrictedForNumberAlphaAndBetaElectrons(int nAlpha, int nBeta) const {
  assert(nAlpha >= 0 && nBeta >= 0);
  assert(coefficientMatrix_.isUnrestricted());
  const auto& cA = coefficientMatrix_.alphaMatrix();
  const auto& cB = coefficientMatrix_.betaMatrix();

  Eigen::MatrixXd alphaMatrix = calculateDensityMatrix(cA, nAlpha);
  Eigen::MatrixXd betaMatrix = calculateDensityMatrix(cB, nBeta);

  DensityMatrix densityMatrix;
  densityMatrix.setDensity(std::move(alphaMatrix), std::move(betaMatrix), nAlpha, nBeta);
  return densityMatrix;
}

DensityMatrix DensityMatrixBuilder::generateRestrictedForSpecifiedOrbitals(const std::vector<int>& doublyOccupiedOrbitals) const {
  assert(coefficientMatrix_.isRestricted());
  const auto& C = coefficientMatrix_.restrictedMatrix();
  auto nAOs = C.rows();
  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(nAOs, nAOs);

  for (auto o : doublyOccupiedOrbitals) {
    assert(o < nAOs && "Orbital index larger than the number of atomic orbitals.");
    P += 2 * calculateSingleOrbitalDensity(C.col(o));
  }

  DensityMatrix densityMatrix;
  auto nSpecifiedOrbitals = static_cast<int>(2 * doublyOccupiedOrbitals.size());
  densityMatrix.setDensity(std::move(P), nSpecifiedOrbitals);
  return densityMatrix;
}

DensityMatrix DensityMatrixBuilder::generateUnrestrictedForSpecifiedOrbitals(const std::vector<int>& alphaOrbitals,
                                                                             const std::vector<int>& betaOrbitals) const {
  assert(coefficientMatrix_.isUnrestricted());
  const auto& Ca = coefficientMatrix_.alphaMatrix();
  const auto& Cb = coefficientMatrix_.betaMatrix();
  auto nAOs = Ca.rows();
  Eigen::MatrixXd Pa = Eigen::MatrixXd::Zero(nAOs, nAOs);
  Eigen::MatrixXd Pb = Eigen::MatrixXd::Zero(nAOs, nAOs);

  for (auto o : alphaOrbitals) {
    assert(o < nAOs && "Orbital index larger than the number of atomic orbitals.");
    Pa += calculateSingleOrbitalDensity(Ca.col(o));
  }
  for (auto o : betaOrbitals) {
    assert(o < nAOs && "Orbital index larger than the number of atomic orbitals.");
    Pb += calculateSingleOrbitalDensity(Cb.col(o));
  }

  DensityMatrix densityMatrix;
  auto nAlpha = static_cast<int>(alphaOrbitals.size());
  auto nBeta = static_cast<int>(betaOrbitals.size());
  densityMatrix.setDensity(std::move(Pa), std::move(Pb), nAlpha, nBeta);
  return densityMatrix;
}

DensityMatrix
DensityMatrixBuilder::generateRestrictedWithSwaps(const std::vector<MolecularOrbitalsManipulation::DeprecatedSwap>& swaps,
                                                  int nElectrons) const {
  assert(coefficientMatrix_.isRestricted());
  assert(nElectrons % 2 == 0 && "The number of electrons is not even.");
  const auto& C = coefficientMatrix_.restrictedMatrix();
  Eigen::MatrixXd P = 2 * calculateDensityMatrix(C, nElectrons / 2);
  P += 2 * calculateDifferenceSwapDensity(C, swaps, nElectrons / 2 - 1);

  DensityMatrix densityMatrix;
  densityMatrix.setDensity(std::move(P), nElectrons);
  return densityMatrix;
}

DensityMatrix DensityMatrixBuilder::generateUnrestrictedWithSwaps(
    const std::vector<MolecularOrbitalsManipulation::DeprecatedSwap>& alphaSwaps,
    const std::vector<MolecularOrbitalsManipulation::DeprecatedSwap>& betaSwaps, int nAlpha, int nBeta) const {
  assert(coefficientMatrix_.isUnrestricted());
  const auto& cA = coefficientMatrix_.alphaMatrix();
  const auto& cB = coefficientMatrix_.betaMatrix();

  Eigen::MatrixXd alphaMatrix = calculateDensityMatrix(cA, nAlpha);
  Eigen::MatrixXd betaMatrix = calculateDensityMatrix(cB, nBeta);

  alphaMatrix += calculateDifferenceSwapDensity(cA, alphaSwaps, nAlpha - 1);
  betaMatrix += calculateDifferenceSwapDensity(cB, betaSwaps, nBeta - 1);

  DensityMatrix densityMatrix;
  densityMatrix.setDensity(std::move(alphaMatrix), std::move(betaMatrix), nAlpha, nBeta);
  return densityMatrix;
}

Eigen::MatrixXd DensityMatrixBuilder::calculateSingleOrbitalDensity(const Eigen::VectorXd& eigenvector) {
  return eigenvector * eigenvector.transpose();
}

Eigen::MatrixXd DensityMatrixBuilder::calculateBlockOrbitalDensity(const Eigen::MatrixXd& eigenvectors) {
  return eigenvectors * eigenvectors.transpose();
}

Eigen::MatrixXd
DensityMatrixBuilder::calculateDifferenceSwapDensity(const Eigen::MatrixXd& coefficientMatrix,
                                                     const std::vector<MolecularOrbitalsManipulation::DeprecatedSwap>& swaps,
                                                     int homoIndex) {
  auto dim = coefficientMatrix.rows();
  Eigen::MatrixXd dP = Eigen::MatrixXd::Zero(dim, dim);

  for (auto s : swaps) {
    auto oldOrbital = homoIndex + 1 + s.lowLyingOrbitalNotToFill_;
    auto newOrbital = homoIndex + s.highLyingOrbitalToFill_;
    dP -= calculateSingleOrbitalDensity(coefficientMatrix.col(oldOrbital));
    dP += calculateSingleOrbitalDensity(coefficientMatrix.col(newOrbital));
  }

  return dP;
}

DensityMatrix DensityMatrixBuilder::generateRestrictedForSpecifiedPartlyOccupiedOrbitals(
    const std::vector<DensityMatrixBuilder::PartlyOccupiedOrbital>& orbitals) {
  assert(coefficientMatrix_.isRestricted());
  const auto& C = coefficientMatrix_.restrictedMatrix();
  auto nAOs = C.rows();
  DensityMatrix densityMatrix;
  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(nAOs, nAOs);
  densityMatrix.setDensity(std::move(P), 0);

  for (auto o : orbitals) {
    assert(o.first < nAOs && "Orbital index larger than the number of atomic orbitals.");
    densityMatrix += generateRestrictedForSpecifiedOrbitals({o.first}) * (o.second / 2);
  }

  return densityMatrix;
}

DensityMatrix DensityMatrixBuilder::generateUnrestrictedForSpecifiedPartlyOccupiedOrbitals(
    const std::vector<DensityMatrixBuilder::PartlyOccupiedOrbital>& alphaOrbitals,
    const std::vector<DensityMatrixBuilder::PartlyOccupiedOrbital>& betaOrbitals) {
  assert(coefficientMatrix_.isUnrestricted());
  const auto& Ca = coefficientMatrix_.alphaMatrix();
  auto nAOs = Ca.rows();
  DensityMatrix densityMatrix;
  Eigen::MatrixXd Pa = Eigen::MatrixXd::Zero(nAOs, nAOs);
  Eigen::MatrixXd Pb = Eigen::MatrixXd::Zero(nAOs, nAOs);
  densityMatrix.setDensity(std::move(Pa), std::move(Pb), 0, 0);

  for (auto o : alphaOrbitals) {
    assert(o.first < nAOs && "Orbital index larger than the number of atomic orbitals.");
    densityMatrix += generateUnrestrictedForSpecifiedOrbitals({o.first}, {}) * o.second;
  }
  for (auto o : betaOrbitals) {
    assert(o.first < nAOs && "Orbital index larger than the number of atomic orbitals.");
    densityMatrix += generateUnrestrictedForSpecifiedOrbitals({}, {o.first}) * o.second;
  }

  return densityMatrix;
}

DensityMatrix
DensityMatrixBuilder::generateRestrictedWithMixing(const std::vector<MolecularOrbitalsManipulation::DeprecatedMix>& mix,
                                                   int nElectrons) const {
  assert(coefficientMatrix_.isRestricted());

  const auto& C = coefficientMatrix_.restrictedMatrix();
  Eigen::MatrixXd P = 2 * calculateDensityMatrix(C, nElectrons / 2);

  P += 2 * calculateDifferenceMixDensity(C, mix, nElectrons / 2 - 1);

  DensityMatrix densityMatrix;
  densityMatrix.setDensity(std::move(P), nElectrons);
  return densityMatrix;
}

DensityMatrix DensityMatrixBuilder::generateUnrestrictedWithMixing(
    const std::vector<MolecularOrbitalsManipulation::DeprecatedMix>& alphaMix,
    const std::vector<MolecularOrbitalsManipulation::DeprecatedMix>& betaMix, int nAlpha, int nBeta) const {
  assert(coefficientMatrix_.isUnrestricted());
  const auto& cA = coefficientMatrix_.alphaMatrix();
  const auto& cB = coefficientMatrix_.betaMatrix();

  Eigen::MatrixXd alphaMatrix = calculateDensityMatrix(cA, nAlpha);
  Eigen::MatrixXd betaMatrix = calculateDensityMatrix(cB, nBeta);

  alphaMatrix += calculateDifferenceMixDensity(cA, alphaMix, nAlpha - 1);
  betaMatrix += calculateDifferenceMixDensity(cB, betaMix, nBeta - 1);

  DensityMatrix densityMatrix;
  densityMatrix.setDensity(std::move(alphaMatrix), std::move(betaMatrix), nAlpha, nBeta);
  return densityMatrix;
}

Eigen::MatrixXd
DensityMatrixBuilder::calculateDifferenceMixDensity(const Eigen::MatrixXd& coefficientMatrix,
                                                    const std::vector<MolecularOrbitalsManipulation::DeprecatedMix>& mix,
                                                    int homoIndex) {
  auto dim = coefficientMatrix.rows();
  Eigen::MatrixXd dP = Eigen::MatrixXd::Zero(dim, dim);

  for (auto m : mix) {
    auto oldOrbital = homoIndex + 1 + m.lowLyingOrbitalNotToFill_;
    auto newOrbital = homoIndex + m.highLyingOrbitalToFill_;
    dP -= calculateSingleOrbitalDensity(coefficientMatrix.col(oldOrbital));
    dP += calculateBlockOrbitalDensity(coefficientMatrix.col(oldOrbital) * std::cos(m.angleInRad_) +
                                       coefficientMatrix.col(newOrbital) * std::sin(m.angleInRad_));
  }

  return dP;
}

} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine
