/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Utils {

using namespace testing;

class AMatrixWithDerivativesTest : public Test {
 public:
  Eigen::MatrixXd randomValuesA, randomValuesB;
  MatrixWithDerivatives matrixA, matrixB;
  Eigen::Matrix<AutomaticDifferentiation::First3D, Eigen::Dynamic, Eigen::Dynamic> randomDerivativesA, randomDerivativesB;
  int arbitraryDimension = 10;

 private:
  void SetUp() final {
    randomValuesA = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    randomValuesB = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    matrixA.reset(arbitraryDimension);
    matrixB.reset(arbitraryDimension);
    randomDerivativesA.resize(arbitraryDimension, arbitraryDimension);
    randomDerivativesB.resize(arbitraryDimension, arbitraryDimension);
    for (int i = 0; i < arbitraryDimension; ++i) {
      for (int j = 0; j < arbitraryDimension; ++j) {
        randomDerivativesA(i, j) = {randomValuesA(i, j), Eigen::Vector3d::Random()};
        randomDerivativesB(i, j) = {randomValuesB(i, j), Eigen::Vector3d::Random()};
      }
    }
    matrixA.get<DerivativeOrder::Zero>() = randomValuesA;
    matrixA.get<DerivativeOrder::One>() = randomDerivativesA;
    matrixB.get<DerivativeOrder::Zero>() = randomValuesB;
    matrixB.get<DerivativeOrder::One>() = randomDerivativesB;
  }
};

TEST_F(AMatrixWithDerivativesTest, CanCorrectlyBeConstructed) {
  MatrixWithDerivatives matrix;
}

TEST_F(AMatrixWithDerivativesTest, CanSumTwoMatricesWithDerivatives) {
  auto result = matrixA + matrixB;
  Eigen::MatrixXd actualResult = randomValuesA + randomValuesB;

  for (int i = 0; i < arbitraryDimension; ++i) {
    for (int j = 0; j < arbitraryDimension; ++j) {
      EXPECT_DOUBLE_EQ(result.get<DerivativeOrder::Zero>()(i, j), actualResult(i, j));
    }
  }
}

TEST_F(AMatrixWithDerivativesTest, CanSubtractTwoMatricesWithDerivatives) {
  auto result = matrixA - matrixB;
  Eigen::MatrixXd actualResult = randomValuesA - randomValuesB;
  for (int i = 0; i < arbitraryDimension; ++i) {
    for (int j = 0; j < arbitraryDimension; ++j) {
      EXPECT_DOUBLE_EQ(result.get<DerivativeOrder::Zero>()(i, j), actualResult(i, j));
    }
  }
}

TEST_F(AMatrixWithDerivativesTest, CanSumDerivativesOfTwoMatricesWithDerivatives) {
  auto result = matrixA + matrixB;
  Eigen::Matrix<AutomaticDifferentiation::First3D, Eigen::Dynamic, Eigen::Dynamic> actualResult =
      randomDerivativesA + randomDerivativesB;

  for (int i = 0; i < arbitraryDimension; ++i) {
    for (int j = 0; j < arbitraryDimension; ++j) {
      EXPECT_DOUBLE_EQ(result.get<DerivativeOrder::One>()(i, j).value(),
                       (randomDerivativesA + randomDerivativesB)(i, j).value());
      EXPECT_TRUE(result.get<DerivativeOrder::One>()(i, j).derivatives().isApprox(
          (randomDerivativesA + randomDerivativesB)(i, j).derivatives()));
    }
  }
}

TEST_F(AMatrixWithDerivativesTest, CanSubtractDerivativesOfTwoMatricesWithDerivatives) {
  auto result = matrixA - matrixB;
  Eigen::Matrix<AutomaticDifferentiation::First3D, Eigen::Dynamic, Eigen::Dynamic> actualResult =
      randomDerivativesA - randomDerivativesB;

  for (int i = 0; i < arbitraryDimension; ++i) {
    for (int j = 0; j < arbitraryDimension; ++j) {
      EXPECT_DOUBLE_EQ(result.get<DerivativeOrder::One>()(i, j).value(), actualResult(i, j).value());
      EXPECT_TRUE(result.get<DerivativeOrder::One>()(i, j).derivatives().isApprox(
          (randomDerivativesA - randomDerivativesB)(i, j).derivatives()));
    }
  }
}

} // namespace Utils
} // namespace Scine
