/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/Scf/LcaoUtils/MolecularOrbitalsManipulation.h>

namespace Scine {
namespace Utils {

namespace LcaoUtils {

MolecularOrbitals
MolecularOrbitalsManipulation::createRestrictedWithSwaps(const MolecularOrbitals& orbitals,
                                                         const std::vector<MolecularOrbitalsManipulation::Swap>& swaps) {
  MolecularOrbitals mo = orbitals;
  swapRestricted(mo, swaps);
  return mo;
}

MolecularOrbitals MolecularOrbitalsManipulation::createUnrestrictedWithSwaps(
    const MolecularOrbitals& orbitals, const std::vector<MolecularOrbitalsManipulation::Swap>& alphaSwaps,
    const std::vector<MolecularOrbitalsManipulation::Swap>& betaSwaps) {
  MolecularOrbitals mo = orbitals;
  swapUnrestricted(mo, alphaSwaps, betaSwaps);
  return mo;
}

MolecularOrbitals
MolecularOrbitalsManipulation::createRestrictedWithMixes(const MolecularOrbitals& orbitals,
                                                         const std::vector<MolecularOrbitalsManipulation::Mix>& mixes) {
  MolecularOrbitals mo = orbitals;
  mixRestricted(mo, mixes);
  return mo;
}

MolecularOrbitals MolecularOrbitalsManipulation::createUnrestrictedWithMixes(
    const MolecularOrbitals& orbitals, const std::vector<MolecularOrbitalsManipulation::Mix>& alphaMixes,
    const std::vector<MolecularOrbitalsManipulation::Mix>& betaMixes) {
  MolecularOrbitals mo = orbitals;
  mixUnrestricted(mo, alphaMixes, betaMixes);
  return mo;
}

void MolecularOrbitalsManipulation::swapRestricted(MolecularOrbitals& orbitals,
                                                   const std::vector<MolecularOrbitalsManipulation::Swap>& swaps) {
  swapOrbitalCoefficients(orbitals.restrictedMatrix(), swaps);
}

void MolecularOrbitalsManipulation::swapUnrestricted(MolecularOrbitals& orbitals,
                                                     const std::vector<MolecularOrbitalsManipulation::Swap>& alphaSwaps,
                                                     const std::vector<MolecularOrbitalsManipulation::Swap>& betaSwaps) {
  swapOrbitalCoefficients(orbitals.alphaMatrix(), alphaSwaps);
  swapOrbitalCoefficients(orbitals.betaMatrix(), betaSwaps);
}

void MolecularOrbitalsManipulation::mixRestricted(MolecularOrbitals& orbitals,
                                                  const std::vector<MolecularOrbitalsManipulation::Mix>& mixes) {
  mixOrbitalCoefficients(orbitals.restrictedMatrix(), mixes);
}

void MolecularOrbitalsManipulation::mixUnrestricted(MolecularOrbitals& orbitals,
                                                    const std::vector<MolecularOrbitalsManipulation::Mix>& alphaMixes,
                                                    const std::vector<MolecularOrbitalsManipulation::Mix>& betaMixes) {
  mixOrbitalCoefficients(orbitals.alphaMatrix(), alphaMixes);
  mixOrbitalCoefficients(orbitals.betaMatrix(), betaMixes);
}

void MolecularOrbitalsManipulation::swapOrbitalCoefficients(Eigen::MatrixXd& coefficientMatrix,
                                                            const std::vector<MolecularOrbitalsManipulation::Swap>& swaps) {
  for (const auto& s : swaps)
    coefficientMatrix.col(s.orbital1_).swap(coefficientMatrix.col(s.orbital2_));
}

void MolecularOrbitalsManipulation::mixOrbitalCoefficients(Eigen::MatrixXd& coefficientMatrix,
                                                           const std::vector<MolecularOrbitalsManipulation::Mix>& mixes) {
  Eigen::VectorXd temporaryOrbitalCoefficients(coefficientMatrix.rows());
  for (const auto& m : mixes) {
    // NB: for the second orbital, we substract the second term to ensure that the two new orbitals will still be
    // orthonormal.
    temporaryOrbitalCoefficients = coefficientMatrix.col(m.orbital1_);
    coefficientMatrix.col(m.orbital1_) = coefficientMatrix.col(m.orbital1_) * std::cos(m.angleInRad_) +
                                         coefficientMatrix.col(m.orbital2_) * std::sin(m.angleInRad_);
    coefficientMatrix.col(m.orbital2_) = coefficientMatrix.col(m.orbital2_) * std::cos(m.angleInRad_) -
                                         temporaryOrbitalCoefficients * std::sin(m.angleInRad_);
  }
}
} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine
