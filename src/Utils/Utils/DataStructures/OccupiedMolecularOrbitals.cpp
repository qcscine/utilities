/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "OccupiedMolecularOrbitals.h"
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>

namespace Scine {
namespace Utils {

OccupiedMolecularOrbitals::OccupiedMolecularOrbitals(const MolecularOrbitals& allOrbitals,
                                                     const LcaoUtils::ElectronicOccupation& occupation) {
  if (allOrbitals.isRestricted()) {
    constructRestricted(allOrbitals, occupation);
  }
  else {
    constructUnrestricted(allOrbitals, occupation);
  }
}

void OccupiedMolecularOrbitals::constructRestricted(const MolecularOrbitals& allOrbitals,
                                                    const LcaoUtils::ElectronicOccupation& occupation) {
  unrestricted_ = false;
  const auto& allOrbitalsMatrix = allOrbitals.restrictedMatrix();
  const auto& filledOrbitals = occupation.getFilledRestrictedOrbitals();

  auto restrictedMatrix = calculateMatrixForFilledOrbitals(allOrbitalsMatrix, filledOrbitals);
  matrix_.setRestrictedMatrix(std::move(restrictedMatrix));
}

void OccupiedMolecularOrbitals::constructUnrestricted(const MolecularOrbitals& allOrbitals,
                                                      const LcaoUtils::ElectronicOccupation& occupation) {
  unrestricted_ = true;
  const auto& fullAlphaMatrix = allOrbitals.alphaMatrix();
  const auto& fullBetaMatrix = allOrbitals.betaMatrix();
  const auto& filledAlphaOrbitals = occupation.getFilledAlphaOrbitals();
  const auto& filledBetaOrbitals = occupation.getFilledBetaOrbitals();

  auto alphaMatrix = calculateMatrixForFilledOrbitals(fullAlphaMatrix, filledAlphaOrbitals);
  auto betaMatrix = calculateMatrixForFilledOrbitals(fullBetaMatrix, filledBetaOrbitals);
  matrix_.setAlphaMatrix(std::move(alphaMatrix));
  matrix_.setBetaMatrix(std::move(betaMatrix));
}

Eigen::MatrixXd OccupiedMolecularOrbitals::calculateMatrixForFilledOrbitals(const Eigen::MatrixXd& matrixWithAllOrbitals,
                                                                            const std::vector<int>& filledOrbitalsIndexes) {
  auto nBasisFunctions = matrixWithAllOrbitals.rows();
  auto nOrbitals = filledOrbitalsIndexes.size();
  Eigen::MatrixXd matrixForFilledOrbitals(nBasisFunctions, nOrbitals);

  int currentColumn = 0;
  for (auto o : filledOrbitalsIndexes) {
    matrixForFilledOrbitals.col(currentColumn) = matrixWithAllOrbitals.col(o);
    currentColumn++;
  }
  return matrixForFilledOrbitals;
}

void OccupiedMolecularOrbitals::makeUnrestricted() {
  if (unrestricted_) {
    return;
  }
  // Set alpha by copy, set beta by move.
  matrix_.setAlphaMatrix(matrix_.restrictedMatrix());
  matrix_.setBetaMatrix(std::move(matrix_.restrictedMatrix()));
  unrestricted_ = true;
}

OccupiedMolecularOrbitals OccupiedMolecularOrbitals::toUnrestricted() const {
  OccupiedMolecularOrbitals mo = *this;
  mo.makeUnrestricted();
  return mo;
}

} // namespace Utils
} // namespace Scine
