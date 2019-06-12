/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

// Alain 15.08.2015: workaround after VS2013 & intel XE 2015
#ifndef MKL_BLAS
#  define MKL_BLAS MKL_DOMAIN_BLAS
#endif

#include <Utils/MethodEssentials/util/LcaoUtil/ElectronicOccupation.h>
#include <Utils/MethodEssentials/util/LcaoUtil/HFWaveFunctionOverlap.h>
#include <Utils/MethodEssentials/util/MolecularOrbitals.h>
#include <Utils/MethodEssentials/util/OccupiedMolecularOrbitals.h>
#include <gmock/gmock.h>
#include <Eigen/Eigenvalues>
#include <cmath>

using namespace testing;
namespace Scine {
namespace Utils {
namespace LcaoUtil {
namespace Tests {

class AHFWaveFunctionOverlap : public Test {
 protected:
  void SetUp() override;

 public:
  MolecularOrbitals randomOrthonormalRhfOrbitals;
  MolecularOrbitals randomNonOrthonormalRhfOrbitals;
  MolecularOrbitals randomOrthonormalUhfOrbitals;
  MolecularOrbitals randomNonOrthonormalUhfOrbitals;
  Eigen::MatrixXd randomOverlapMatrixForRhf;
  Eigen::MatrixXd randomOverlapMatrixForUhf;
  const int nOrbitals = 8;

  Eigen::MatrixXd generateOrthogonalCoefficientMatrix();
  Eigen::MatrixXd generateNonOrthogonalCoefficientMatrix(const Eigen::MatrixXd& overlap);
  Eigen::MatrixXd generateRandomOverlapMatrix();
};

void AHFWaveFunctionOverlap::SetUp() {
  randomOverlapMatrixForRhf = generateRandomOverlapMatrix();
  randomOverlapMatrixForUhf = generateRandomOverlapMatrix();
  auto randomOrthogonalCoefficientMatrix1 = generateOrthogonalCoefficientMatrix();
  auto randomOrthogonalCoefficientMatrix2 = generateOrthogonalCoefficientMatrix();
  auto randomOrthogonalCoefficientMatrix3 = generateOrthogonalCoefficientMatrix();
  auto randomNonOrthogonalCoefficientMatrix1 = generateNonOrthogonalCoefficientMatrix(randomOverlapMatrixForRhf);
  auto randomNonOrthogonalCoefficientMatrix2 = generateNonOrthogonalCoefficientMatrix(randomOverlapMatrixForUhf);
  auto randomNonOrthogonalCoefficientMatrix3 = generateNonOrthogonalCoefficientMatrix(randomOverlapMatrixForUhf);
  randomOrthonormalRhfOrbitals = MolecularOrbitals::createFromRestrictedCoefficients(randomOrthogonalCoefficientMatrix1);
  randomNonOrthonormalRhfOrbitals = MolecularOrbitals::createFromRestrictedCoefficients(randomNonOrthogonalCoefficientMatrix1);
  randomOrthonormalUhfOrbitals = MolecularOrbitals::createFromUnrestrictedCoefficients(randomOrthogonalCoefficientMatrix2,
                                                                                       randomOrthogonalCoefficientMatrix3);
  randomNonOrthonormalUhfOrbitals = MolecularOrbitals::createFromUnrestrictedCoefficients(
      randomNonOrthogonalCoefficientMatrix2, randomNonOrthogonalCoefficientMatrix3);
}

Eigen::MatrixXd AHFWaveFunctionOverlap::generateRandomOverlapMatrix() {
  Eigen::MatrixXd random = Eigen::MatrixXd::Random(nOrbitals, nOrbitals);
  random /= 10;
  for (int i = 0; i < nOrbitals; ++i)
    random(i, i) = 1;
  return .5 * (random + random.transpose());
}

Eigen::MatrixXd AHFWaveFunctionOverlap::generateOrthogonalCoefficientMatrix() {
  Eigen::MatrixXd matrixToDiagonalize = Eigen::MatrixXd::Random(nOrbitals, nOrbitals);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  es.compute(matrixToDiagonalize);
  return es.eigenvectors();
}

Eigen::MatrixXd AHFWaveFunctionOverlap::generateNonOrthogonalCoefficientMatrix(const Eigen::MatrixXd& overlap) {
  Eigen::MatrixXd matrixToDiagonalize = 20 * Eigen::MatrixXd::Random(nOrbitals, nOrbitals);
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es;
  es.compute(matrixToDiagonalize, overlap);
  return es.eigenvectors();
}

TEST_F(AHFWaveFunctionOverlap, returnsOneForIdenticalRhfWaveFunctions) {
  LcaoUtil::ElectronicOccupation rhfOccupation;
  rhfOccupation.fillSpecifiedRestrictedOrbitals({0, 2, 3});
  OccupiedMolecularOrbitals orthonormalRhfOrbitals{randomOrthonormalRhfOrbitals, rhfOccupation};
  OccupiedMolecularOrbitals nonOrthonormalRhfOrbitals{randomNonOrthonormalRhfOrbitals, rhfOccupation};

  double selfRhfOrthonormal =
      HFWaveFunctionOverlap::calculateOrthonormalOverlap(orthonormalRhfOrbitals, orthonormalRhfOrbitals);
  double selfRhfNonOrthonormal = HFWaveFunctionOverlap::calculateNonOrthonormalOverlap(
      nonOrthonormalRhfOrbitals, nonOrthonormalRhfOrbitals, randomOverlapMatrixForRhf);
  ASSERT_THAT(selfRhfOrthonormal, DoubleNear(1, 1e-12));
  ASSERT_THAT(selfRhfNonOrthonormal, DoubleNear(1, 1e-12));
}

TEST_F(AHFWaveFunctionOverlap, returnsOneForIdenticalUhfWaveFunctions) {
  LcaoUtil::ElectronicOccupation uhfOccupation;
  uhfOccupation.fillSpecifiedUnrestrictedOrbitals({1, 2, 4}, {5, 6});
  OccupiedMolecularOrbitals orthonormalUhfOrbitals{randomOrthonormalUhfOrbitals, uhfOccupation};
  OccupiedMolecularOrbitals nonOrthonormalUhfOrbitals{randomNonOrthonormalUhfOrbitals, uhfOccupation};

  double selfUhfOrthonormal =
      HFWaveFunctionOverlap::calculateOrthonormalOverlap(orthonormalUhfOrbitals, orthonormalUhfOrbitals);
  double selfUhfNonOrthonormal = HFWaveFunctionOverlap::calculateNonOrthonormalOverlap(
      nonOrthonormalUhfOrbitals, nonOrthonormalUhfOrbitals, randomOverlapMatrixForUhf);
  ASSERT_THAT(selfUhfOrthonormal, DoubleNear(1, 1e-12));
  ASSERT_THAT(selfUhfNonOrthonormal, DoubleNear(1, 1e-12));
}

TEST_F(AHFWaveFunctionOverlap, returnsZeroIfOneOrbitalIsSwitched) {
  LcaoUtil::ElectronicOccupation uhfOccupation1;
  LcaoUtil::ElectronicOccupation uhfOccupation2;
  uhfOccupation1.fillSpecifiedUnrestrictedOrbitals({1, 2, 4}, {5, 6});
  uhfOccupation2.fillSpecifiedUnrestrictedOrbitals({1, 2, 4}, {5, 7});
  OccupiedMolecularOrbitals orthonormalUhfOrbitals1{randomOrthonormalUhfOrbitals, uhfOccupation1};
  OccupiedMolecularOrbitals nonOrthonormalUhfOrbitals1{randomNonOrthonormalUhfOrbitals, uhfOccupation1};
  OccupiedMolecularOrbitals orthonormalUhfOrbitals2{randomOrthonormalUhfOrbitals, uhfOccupation2};
  OccupiedMolecularOrbitals nonOrthonormalUhfOrbitals2{randomNonOrthonormalUhfOrbitals, uhfOccupation2};

  double selfUhfOrthonormal =
      HFWaveFunctionOverlap::calculateOrthonormalOverlap(orthonormalUhfOrbitals1, orthonormalUhfOrbitals2);
  double selfUhfNonOrthonormal = HFWaveFunctionOverlap::calculateNonOrthonormalOverlap(
      nonOrthonormalUhfOrbitals1, nonOrthonormalUhfOrbitals2, randomOverlapMatrixForUhf);
  ASSERT_THAT(selfUhfOrthonormal, DoubleNear(0, 1e-12));
  ASSERT_THAT(selfUhfNonOrthonormal, DoubleNear(0, 1e-12));
}

TEST_F(AHFWaveFunctionOverlap, ReturnsOneForOverlapBetweenRhfAndCorrespondingUhfWaveFunctions) {
  LcaoUtil::ElectronicOccupation rhfOccupation;
  LcaoUtil::ElectronicOccupation uhfOccupation;
  rhfOccupation.fillSpecifiedRestrictedOrbitals({0, 2, 3});
  uhfOccupation.fillSpecifiedUnrestrictedOrbitals({0, 2, 3}, {0, 2, 3});

  OccupiedMolecularOrbitals orthonormalRhfOrbitals{randomOrthonormalRhfOrbitals, rhfOccupation};
  OccupiedMolecularOrbitals nonOrthonormalRhfOrbitals{randomNonOrthonormalRhfOrbitals, rhfOccupation};
  auto orthonormalUhfOrbitals = orthonormalRhfOrbitals.toUnrestricted();
  auto nonOrthonormalUhfOrbitals = nonOrthonormalRhfOrbitals.toUnrestricted();

  double orthonormalDistance =
      HFWaveFunctionOverlap::calculateOrthonormalOverlap(orthonormalRhfOrbitals, orthonormalUhfOrbitals);
  double nonOrthonormalDistance = HFWaveFunctionOverlap::calculateNonOrthonormalOverlap(
      nonOrthonormalRhfOrbitals, nonOrthonormalUhfOrbitals, randomOverlapMatrixForRhf);
  ASSERT_THAT(orthonormalDistance, DoubleNear(1, 1e-12));
  ASSERT_THAT(nonOrthonormalDistance, DoubleNear(1, 1e-12));
}

TEST_F(AHFWaveFunctionOverlap, returnsZeroIfDifferentNumberOfElectrons) {
  LcaoUtil::ElectronicOccupation uhfOccupation1;
  LcaoUtil::ElectronicOccupation uhfOccupation2;
  uhfOccupation1.fillSpecifiedUnrestrictedOrbitals({1, 2, 4}, {1, 2});
  uhfOccupation2.fillSpecifiedUnrestrictedOrbitals({1, 2}, {1, 2, 4});
  OccupiedMolecularOrbitals orthonormalUhfOrbitals1{randomOrthonormalUhfOrbitals, uhfOccupation1};
  OccupiedMolecularOrbitals nonOrthonormalUhfOrbitals1{randomNonOrthonormalUhfOrbitals, uhfOccupation1};
  OccupiedMolecularOrbitals orthonormalUhfOrbitals2{randomOrthonormalUhfOrbitals, uhfOccupation2};
  OccupiedMolecularOrbitals nonOrthonormalUhfOrbitals2{randomNonOrthonormalUhfOrbitals, uhfOccupation2};

  ASSERT_ANY_THROW(HFWaveFunctionOverlap::calculateOrthonormalOverlap(orthonormalUhfOrbitals1, orthonormalUhfOrbitals2));
  ASSERT_ANY_THROW(HFWaveFunctionOverlap::calculateNonOrthonormalOverlap(
      nonOrthonormalUhfOrbitals1, nonOrthonormalUhfOrbitals2, randomOverlapMatrixForUhf));
}
} // namespace Tests
} // namespace LcaoUtil
} // namespace Utils
} // namespace Scine
