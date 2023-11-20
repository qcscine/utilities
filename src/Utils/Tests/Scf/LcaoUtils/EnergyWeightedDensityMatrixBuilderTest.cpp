/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "gmock/gmock.h"
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <Utils/Scf/LcaoUtils/EnergyWeightedDensityMatrixBuilder.h>
#include <cmath>

using namespace testing;
namespace Scine {
namespace Utils {
namespace LcaoUtils {
namespace Tests {

class AEnergyWeightedDensityMatrixBuilderTest : public Test {
 public:
  AEnergyWeightedDensityMatrixBuilderTest();

  int dim = 4;
  MolecularOrbitals randomRHFMatrix;
  MolecularOrbitals randomUHFMatrix;
  SingleParticleEnergies randomRHFEnergies;
  SingleParticleEnergies randomUHFEnergies;

  void assertMatrixContainsOnlyZeros(const Eigen::MatrixXd& d);
  void assertMatricesIdentical(const Eigen::MatrixXd& d1, const Eigen::MatrixXd& d2);
  void testReturnsZeroForRHF();
  void testReturnsZeroForUHF();
  void testOrbitalsCanBeSpecifiedForRHF();
  void testOrbitalsCanBeSpecifiedForUHF();
};

AEnergyWeightedDensityMatrixBuilderTest::AEnergyWeightedDensityMatrixBuilderTest() {
  Eigen::MatrixXd m1 = Eigen::MatrixXd::Random(dim, dim);
  Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(dim, dim);
  Eigen::MatrixXd m3 = Eigen::MatrixXd::Random(dim, dim);
  randomRHFMatrix = MolecularOrbitals::createFromRestrictedCoefficients(std::move(m1));
  randomUHFMatrix = MolecularOrbitals::createFromUnrestrictedCoefficients(std::move(m2), std::move(m3));

  Eigen::Vector4d v1 = Eigen::Vector4d::Random();
  Eigen::Vector4d v2 = Eigen::Vector4d::Random();
  Eigen::Vector4d v3 = Eigen::Vector4d::Random();
  randomRHFEnergies.setRestricted(v1);
  randomUHFEnergies.setUnrestricted(v2, v3);
}

void AEnergyWeightedDensityMatrixBuilderTest::testReturnsZeroForUHF() {
  EnergyWeightedDensityMatrixBuilder matrixBuilderUHF(randomUHFMatrix, randomUHFEnergies);
  auto dUHF = matrixBuilderUHF.generateUnrestrictedForNumberAlphaAndBetaElectrons(0, 0);
  assertMatrixContainsOnlyZeros(dUHF);
}

void AEnergyWeightedDensityMatrixBuilderTest::testReturnsZeroForRHF() {
  EnergyWeightedDensityMatrixBuilder matrixBuilderRHF(randomRHFMatrix, randomRHFEnergies);
  auto dRHF = matrixBuilderRHF.generateRestrictedForNumberElectrons(0);
  assertMatrixContainsOnlyZeros(dRHF);
}

void AEnergyWeightedDensityMatrixBuilderTest::assertMatrixContainsOnlyZeros(const Eigen::MatrixXd& d) {
  ASSERT_TRUE(d.isZero(1e-10));
}

void AEnergyWeightedDensityMatrixBuilderTest::testOrbitalsCanBeSpecifiedForRHF() {
  EnergyWeightedDensityMatrixBuilder matrixBuilderRHF(randomRHFMatrix, randomRHFEnergies);
  std::vector<int> orbitals = {1, 2};

  auto d = matrixBuilderRHF.generateRestrictedForSpecifiedOrbitals(orbitals);

  const Eigen::MatrixXd& C = randomRHFMatrix.restrictedMatrix();
  Eigen::MatrixXd P = 2 * randomRHFEnergies.getRestrictedLevelEnergy(1) * C.col(1) * C.col(1).transpose();
  P += 2 * randomRHFEnergies.getRestrictedLevelEnergy(2) * C.col(2) * C.col(2).transpose();

  assertMatricesIdentical(d, P);
}

void AEnergyWeightedDensityMatrixBuilderTest::testOrbitalsCanBeSpecifiedForUHF() {
  EnergyWeightedDensityMatrixBuilder matrixBuilderUHF(randomUHFMatrix, randomUHFEnergies);
  std::vector<int> alphaOrbitals = {1, 2};
  std::vector<int> betaOrbitals = {0};

  auto d = matrixBuilderUHF.generateUnrestrictedForSpecifiedOrbitals(alphaOrbitals, betaOrbitals);

  Eigen::MatrixXd Ca = randomUHFMatrix.alphaMatrix();
  Eigen::MatrixXd Cb = randomUHFMatrix.betaMatrix();
  Eigen::MatrixXd P = Ca.col(1) * Ca.col(1).transpose() * randomUHFEnergies.getAlphaLevelEnergy(1);
  P += Ca.col(2) * Ca.col(2).transpose() * randomUHFEnergies.getAlphaLevelEnergy(2);
  P += Cb.col(0) * Cb.col(0).transpose() * randomUHFEnergies.getBetaLevelEnergy(0);

  assertMatricesIdentical(d, P);
}

void AEnergyWeightedDensityMatrixBuilderTest::assertMatricesIdentical(const Eigen::MatrixXd& d1, const Eigen::MatrixXd& d2) {
  ASSERT_TRUE(d1.isApprox(d2, 1e-10));
}

TEST_F(AEnergyWeightedDensityMatrixBuilderTest, ReturnsMatrixFilledWithZerosForIfNoElectrons) {
  testReturnsZeroForRHF();
  testReturnsZeroForUHF();
}

TEST_F(AEnergyWeightedDensityMatrixBuilderTest, AllowsToSpecifyWhichOrbitalsToConsider) {
  testOrbitalsCanBeSpecifiedForRHF();
  testOrbitalsCanBeSpecifiedForUHF();
}

TEST_F(AEnergyWeightedDensityMatrixBuilderTest, SpecifyingLowestOrbitalsIsTheSameAsSpecifyingTheNumberOfElectrons) {
  EnergyWeightedDensityMatrixBuilder matrixBuilderRHF(randomRHFMatrix, randomRHFEnergies);
  EnergyWeightedDensityMatrixBuilder matrixBuilderUHF(randomUHFMatrix, randomUHFEnergies);

  assertMatricesIdentical(matrixBuilderRHF.generateRestrictedForNumberElectrons(6),
                          matrixBuilderRHF.generateRestrictedForSpecifiedOrbitals({0, 1, 2}));
  assertMatricesIdentical(matrixBuilderUHF.generateUnrestrictedForNumberAlphaAndBetaElectrons(4, 3),
                          matrixBuilderUHF.generateUnrestrictedForSpecifiedOrbitals({0, 1, 2, 3}, {0, 1, 2}));
}

} // namespace Tests
} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine
