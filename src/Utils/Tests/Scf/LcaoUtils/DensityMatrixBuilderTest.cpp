/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/Scf/LcaoUtils/DensityMatrixBuilder.h>
#include <gmock/gmock.h>
#include <cmath>

using namespace testing;
namespace Scine {
namespace Utils {
namespace LcaoUtils {
namespace Tests {

class ADensityMatrixBuilderTest : public Test {
 public:
  ADensityMatrixBuilderTest();

  int dim = 4;
  MolecularOrbitals randomRHFMatrix;
  MolecularOrbitals randomUHFMatrix;

  void assertMatrixContainsOnlyZeros(const DensityMatrix& d);
  void assertMatricesIdentical(const DensityMatrix& d1, const DensityMatrix& d2);
  void testReturnsZeroForRHF();
  void testReturnsZeroForUHF();
  void testProducesCorrectNumberOfElectronsForRHF();
  void testProducesCorrectNumberOfElectronsForUHF();
  void testProducesCorrectNumberOfAOsForRHF();
  void testProducesCorrectNumberOfAOsForUHF();
  void testOrbitalsCanBeSpecifiedForRHF();
  void testOrbitalsCanBeSpecifiedForUHF();
  void testSwapsCanBeSpecifiedForRHF();
  void testSwapsCanBeSpecifiedForUHF();
  void testPartialOccupationForRHF();
  void testPartialOccupationForUHF();
  void testOrbitalMixingForRHF();
  void testOrbitalMixingForUHF();
};

ADensityMatrixBuilderTest::ADensityMatrixBuilderTest() {
  Eigen::MatrixXd m1 = Eigen::MatrixXd::Random(dim, dim);
  Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(dim, dim);
  Eigen::MatrixXd m3 = Eigen::MatrixXd::Random(dim, dim);
  randomRHFMatrix = MolecularOrbitals::createFromRestrictedCoefficients(std::move(m1));
  randomUHFMatrix = MolecularOrbitals::createFromUnrestrictedCoefficients(std::move(m2), std::move(m3));
}

void ADensityMatrixBuilderTest::testReturnsZeroForUHF() {
  DensityMatrixBuilder matrixBuilderUHF(randomUHFMatrix);
  auto dUHF = matrixBuilderUHF.generateUnrestrictedForNumberAlphaAndBetaElectrons(0, 0);
  assertMatrixContainsOnlyZeros(dUHF);
}

void ADensityMatrixBuilderTest::testReturnsZeroForRHF() {
  DensityMatrixBuilder matrixBuilderRHF(randomRHFMatrix);
  auto dRHF = matrixBuilderRHF.generateRestrictedForNumberElectrons(0);
  assertMatrixContainsOnlyZeros(dRHF);
}

void ADensityMatrixBuilderTest::assertMatrixContainsOnlyZeros(const DensityMatrix& d) {
  if (d.unrestricted()) {
    ASSERT_TRUE(d.alphaMatrix().isZero(1e-10));
    ASSERT_TRUE(d.betaMatrix().isZero(1e-10));
  }
  else {
    ASSERT_TRUE(d.restrictedMatrix().isZero(1e-10));
  }
}

void ADensityMatrixBuilderTest::testProducesCorrectNumberOfElectronsForUHF() {
  DensityMatrixBuilder matrixBuilderUHF(randomUHFMatrix);
  int nAlpha = 3, nBeta = 4;

  auto dUHF = matrixBuilderUHF.generateUnrestrictedForNumberAlphaAndBetaElectrons(nAlpha, nBeta);

  ASSERT_THAT(dUHF.numberElectronsInAlphaMatrix(), Eq(nAlpha));
  ASSERT_THAT(dUHF.numberElectronsInBetaMatrix(), Eq(nBeta));
}

void ADensityMatrixBuilderTest::testProducesCorrectNumberOfAOsForRHF() {
  DensityMatrixBuilder matrixBuilderRHF(randomRHFMatrix);
  int AOs = randomRHFMatrix.restrictedMatrix().rows();

  auto dRHF = matrixBuilderRHF.generateRestrictedForNumberElectrons(6);

  ASSERT_THAT(dRHF.restrictedMatrix().rows(), Eq(AOs));
  ASSERT_THAT(dRHF.restrictedMatrix().cols(), Eq(AOs));
}

void ADensityMatrixBuilderTest::testProducesCorrectNumberOfAOsForUHF() {
  DensityMatrixBuilder matrixBuilderUHF(randomUHFMatrix);
  int AOs = randomUHFMatrix.restrictedMatrix().rows();

  auto dUHF = matrixBuilderUHF.generateUnrestrictedForNumberAlphaAndBetaElectrons(6, 6);

  ASSERT_THAT(dUHF.restrictedMatrix().rows(), Eq(AOs));
  ASSERT_THAT(dUHF.restrictedMatrix().cols(), Eq(AOs));
  ASSERT_THAT(dUHF.alphaMatrix().rows(), Eq(AOs));
  ASSERT_THAT(dUHF.alphaMatrix().cols(), Eq(AOs));
  ASSERT_THAT(dUHF.betaMatrix().rows(), Eq(AOs));
  ASSERT_THAT(dUHF.betaMatrix().cols(), Eq(AOs));
}

void ADensityMatrixBuilderTest::testProducesCorrectNumberOfElectronsForRHF() {
  DensityMatrixBuilder matrixBuilderRHF(randomRHFMatrix);
  int n = 6;

  auto dRHF = matrixBuilderRHF.generateRestrictedForNumberElectrons(n);

  ASSERT_THAT(dRHF.numberElectronsInAlphaMatrix() + dRHF.numberElectronsInBetaMatrix(), Eq(n));
}

void ADensityMatrixBuilderTest::testOrbitalsCanBeSpecifiedForRHF() {
  DensityMatrixBuilder matrixBuilderRHF(randomRHFMatrix);
  std::vector<int> orbitals = {1, 2};

  auto d = matrixBuilderRHF.generateRestrictedForSpecifiedOrbitals(orbitals);

  Eigen::MatrixXd C = randomRHFMatrix.restrictedMatrix();
  Eigen::MatrixXd P = 2 * C.col(1) * C.col(1).transpose();
  P += 2 * C.col(2) * C.col(2).transpose();
  DensityMatrix expected;
  expected.setDensity(std::move(P), 4);

  assertMatricesIdentical(d, expected);
}

void ADensityMatrixBuilderTest::testOrbitalsCanBeSpecifiedForUHF() {
  DensityMatrixBuilder matrixBuilderUHF(randomUHFMatrix);
  std::vector<int> alphaOrbitals = {1, 2};
  std::vector<int> betaOrbitals = {0};

  auto d = matrixBuilderUHF.generateUnrestrictedForSpecifiedOrbitals(alphaOrbitals, betaOrbitals);

  Eigen::MatrixXd Ca = randomUHFMatrix.alphaMatrix();
  Eigen::MatrixXd Cb = randomUHFMatrix.betaMatrix();
  Eigen::MatrixXd Pa = Ca.col(1) * Ca.col(1).transpose();
  Pa += Ca.col(2) * Ca.col(2).transpose();
  Eigen::MatrixXd Pb = Cb.col(0) * Cb.col(0).transpose();
  DensityMatrix expected;
  expected.setDensity(std::move(Pa), std::move(Pb), 2, 1);

  assertMatricesIdentical(d, expected);
}

void ADensityMatrixBuilderTest::assertMatricesIdentical(const DensityMatrix& d1, const DensityMatrix& d2) {
  ASSERT_THAT(d1.getAlphaOccupation(), DoubleNear(d2.getAlphaOccupation(), 1e-12));
  ASSERT_THAT(d1.getBetaOccupation(), DoubleNear(d2.getBetaOccupation(), 1e-12));
  if (d1.unrestricted()) {
    ASSERT_TRUE(d1.betaMatrix().isApprox(d2.betaMatrix(), 1e-10));
    ASSERT_TRUE(d1.alphaMatrix().isApprox(d2.alphaMatrix(), 1e-10));
  }
  else
    ASSERT_TRUE(d1.restrictedMatrix().isApprox(d2.restrictedMatrix(), 1e-10));
}

void ADensityMatrixBuilderTest::testSwapsCanBeSpecifiedForRHF() {
  DensityMatrixBuilder matrixBuilderRHF(randomRHFMatrix);

  std::vector<MolecularOrbitalsManipulation::DeprecatedSwap> swaps = {{-2, 1}, {-1, 2}};
  auto d1 = matrixBuilderRHF.generateRestrictedWithSwaps(swaps, 4);

  auto d2 = matrixBuilderRHF.generateRestrictedForSpecifiedOrbitals({2, 3});

  assertMatricesIdentical(d1, d2);
}

void ADensityMatrixBuilderTest::testSwapsCanBeSpecifiedForUHF() {
  DensityMatrixBuilder matrixBuilderUHF(randomUHFMatrix);

  std::vector<MolecularOrbitalsManipulation::DeprecatedSwap> alphaSwaps = {{-2, 1}, {-1, 2}};
  std::vector<MolecularOrbitalsManipulation::DeprecatedSwap> betaSwaps = {{-1, 1}};
  auto d1 = matrixBuilderUHF.generateUnrestrictedWithSwaps(alphaSwaps, betaSwaps, 2, 1);

  std::vector<int> alphaOrbitals = {2, 3};
  std::vector<int> betaOrbitals = {1};
  auto d2 = matrixBuilderUHF.generateUnrestrictedForSpecifiedOrbitals(alphaOrbitals, betaOrbitals);

  assertMatricesIdentical(d1, d2);
}

void ADensityMatrixBuilderTest::testPartialOccupationForRHF() {
  DensityMatrixBuilder matrixBuilderRHF(randomRHFMatrix);

  double occupation = 0.5876;
  auto d = matrixBuilderRHF.generateRestrictedForSpecifiedPartlyOccupiedOrbitals({{0, 2}, {1, occupation}});

  auto firstPart = matrixBuilderRHF.generateRestrictedForSpecifiedOrbitals({0});
  auto secondPart = matrixBuilderRHF.generateRestrictedForSpecifiedOrbitals({1});
  auto expected = firstPart + secondPart * (occupation / 2.);

  assertMatricesIdentical(d, expected);
}

void ADensityMatrixBuilderTest::testPartialOccupationForUHF() {
  DensityMatrixBuilder matrixBuilderUHF(randomUHFMatrix);

  double occupationA = 0.5876;
  double occupationB = 0.1234;

  auto d = matrixBuilderUHF.generateUnrestrictedForSpecifiedPartlyOccupiedOrbitals({{2, 1}, {3, occupationA}},
                                                                                   {{0, occupationB}});

  auto firstPart = matrixBuilderUHF.generateUnrestrictedForSpecifiedOrbitals({2}, {});
  auto secondPart = matrixBuilderUHF.generateUnrestrictedForSpecifiedOrbitals({3}, {}) * occupationA;
  auto thirdPart = matrixBuilderUHF.generateUnrestrictedForSpecifiedOrbitals({}, {0}) * occupationB;
  auto expected = firstPart + secondPart + thirdPart;

  assertMatricesIdentical(d, expected);
}

void ADensityMatrixBuilderTest::testOrbitalMixingForRHF() {
  DensityMatrixBuilder matrixBuilderRHF(randomRHFMatrix);

  double mixPhase = 0.1234;
  std::vector<MolecularOrbitalsManipulation::DeprecatedMix> mix = {{-2, 1, mixPhase}};
  auto d1 = matrixBuilderRHF.generateRestrictedWithMixing(mix, 4);

  Eigen::MatrixXd mixedRHF = randomRHFMatrix.restrictedMatrix();
  mixedRHF.col(0) = randomRHFMatrix.restrictedMatrix().col(0) * std::cos(mixPhase) +
                    randomRHFMatrix.restrictedMatrix().col(2) * std::sin(mixPhase);
  auto mixedCoefficientRHFMatrix = MolecularOrbitals::createFromRestrictedCoefficients(std::move(mixedRHF));
  DensityMatrixBuilder mixedBuilder(mixedCoefficientRHFMatrix);
  auto expected = mixedBuilder.generateRestrictedForNumberElectrons(4);

  assertMatricesIdentical(d1, expected);
}

void ADensityMatrixBuilderTest::testOrbitalMixingForUHF() {
  DensityMatrixBuilder matrixBuilderUHF(randomUHFMatrix);

  double mixPhaseA = 0.1234;
  double mixPhaseB = 0.8234;
  std::vector<MolecularOrbitalsManipulation::DeprecatedMix> mixA = {{-2, 1, mixPhaseA}};
  std::vector<MolecularOrbitalsManipulation::DeprecatedMix> mixB = {{-1, 1, mixPhaseB}};
  auto d1 = matrixBuilderUHF.generateUnrestrictedWithMixing(mixA, mixB, 2, 1);

  Eigen::MatrixXd mixedAlpha = randomUHFMatrix.alphaMatrix();
  Eigen::MatrixXd mixedBeta = randomUHFMatrix.betaMatrix();
  mixedAlpha.col(0) = randomUHFMatrix.alphaMatrix().col(0) * std::cos(mixPhaseA) +
                      randomUHFMatrix.alphaMatrix().col(2) * std::sin(mixPhaseA);
  mixedBeta.col(0) = randomUHFMatrix.betaMatrix().col(0) * std::cos(mixPhaseB) +
                     randomUHFMatrix.betaMatrix().col(1) * std::sin(mixPhaseB);
  auto mixedCoefficientUHFMatrix = MolecularOrbitals::createFromUnrestrictedCoefficients(mixedAlpha, mixedBeta);
  DensityMatrixBuilder mixedBuilder(mixedCoefficientUHFMatrix);
  auto expected = mixedBuilder.generateUnrestrictedForNumberAlphaAndBetaElectrons(2, 1);

  assertMatricesIdentical(d1, expected);
}

TEST_F(ADensityMatrixBuilderTest, ReturnsMatrixFilledWithZerosForIfNoElectrons) {
  testReturnsZeroForRHF();
  testReturnsZeroForUHF();
}

TEST_F(ADensityMatrixBuilderTest, ProducesADensityMatrixWithCorrectNumberOfElectrons) {
  testProducesCorrectNumberOfElectronsForRHF();
  testProducesCorrectNumberOfElectronsForUHF();
}

TEST_F(ADensityMatrixBuilderTest, AllowsToSpecifyWhichOrbitalsToConsider) {
  testOrbitalsCanBeSpecifiedForRHF();
  testOrbitalsCanBeSpecifiedForUHF();
}

TEST_F(ADensityMatrixBuilderTest, SpecifyingLowestOrbitalsIsTheSameAsSpecifyingTheNumberOfElectrons) {
  DensityMatrixBuilder matrixBuilderRHF(randomRHFMatrix);
  DensityMatrixBuilder matrixBuilderUHF(randomUHFMatrix);

  assertMatricesIdentical(matrixBuilderRHF.generateRestrictedForNumberElectrons(6),
                          matrixBuilderRHF.generateRestrictedForSpecifiedOrbitals({0, 1, 2}));
  assertMatricesIdentical(matrixBuilderUHF.generateUnrestrictedForNumberAlphaAndBetaElectrons(4, 3),
                          matrixBuilderUHF.generateUnrestrictedForSpecifiedOrbitals({0, 1, 2, 3}, {0, 1, 2}));
}

TEST_F(ADensityMatrixBuilderTest, AllowsToSpecifyOrbitalSwaps) {
  testSwapsCanBeSpecifiedForRHF();
  testSwapsCanBeSpecifiedForUHF();
}

TEST_F(ADensityMatrixBuilderTest, AllowsForPartialOccupation) {
  testPartialOccupationForRHF();
  testPartialOccupationForUHF();
}

TEST_F(ADensityMatrixBuilderTest, AllowsForOrbitalMixing) {
  testOrbitalMixingForRHF();
  testOrbitalMixingForUHF();
}

} // namespace Tests
} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine
