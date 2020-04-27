/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/DipoleMatrix.h>
#include <Utils/DataStructures/OctupoleMatrix.h>
#include <Utils/DataStructures/QuadrupoleMatrix.h>
#include <Utils/Math/AutomaticDifferentiation/First3D.h>
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Utils {
namespace Tests {

class AMultipoleMatrixTest : public Test {
 public:
  int arbitraryDimension = 10;
  DipoleMatrix dipoleMatrix;
  QuadrupoleMatrix quadrupoleMatrix;
  OctupoleMatrix octupoleMatrix;
  Eigen::MatrixXd randomXMatrix, randomYMatrix, randomZMatrix, randomXYMatrix, randomXZMatrix, randomYZMatrix,
      randomXYYMatrix, randomXZZMatrix, randomYYZMatrix, randomYZZMatrix;
  Eigen::Matrix<AutomaticDifferentiation::First3D, Eigen::Dynamic, Eigen::Dynamic> randomXderivatives,
      randomYderivatives, randomZderivatives, randomXYderivatives, randomXZderivatives, randomYZderivatives,
      randomXYYderivatives, randomXZZderivatives, randomYYZderivatives, randomYZZderivatives;

 private:
  void SetUp() final {
    randomXMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    randomYMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    randomZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    randomXYMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    randomXZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    randomYZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    randomXYYMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    randomXZZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    randomYYZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    randomYZZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);

    // Generate derivatives
    randomXderivatives.resize(arbitraryDimension, arbitraryDimension);
    randomYderivatives.resize(arbitraryDimension, arbitraryDimension);
    randomZderivatives.resize(arbitraryDimension, arbitraryDimension);
    randomXYderivatives.resize(arbitraryDimension, arbitraryDimension);
    randomXZderivatives.resize(arbitraryDimension, arbitraryDimension);
    randomYZderivatives.resize(arbitraryDimension, arbitraryDimension);
    randomXYYderivatives.resize(arbitraryDimension, arbitraryDimension);
    randomXZZderivatives.resize(arbitraryDimension, arbitraryDimension);
    randomYYZderivatives.resize(arbitraryDimension, arbitraryDimension);
    randomYZZderivatives.resize(arbitraryDimension, arbitraryDimension);
    for (int i = 0; i < arbitraryDimension; ++i) {
      for (int j = 0; j < arbitraryDimension; ++j) {
        randomXderivatives(i, j) = {randomXMatrix(i, j), Eigen::Vector3d::Random()};
        randomYderivatives(i, j) = {randomYMatrix(i, j), Eigen::Vector3d::Random()};
        randomZderivatives(i, j) = {randomZMatrix(i, j), Eigen::Vector3d::Random()};
        randomXYderivatives(i, j) = {randomXYMatrix(i, j), Eigen::Vector3d::Random()};
        randomXZderivatives(i, j) = {randomXZMatrix(i, j), Eigen::Vector3d::Random()};
        randomYZderivatives(i, j) = {randomYZMatrix(i, j), Eigen::Vector3d::Random()};
        randomXYYderivatives(i, j) = {randomXYYMatrix(i, j), Eigen::Vector3d::Random()};
        randomXZZderivatives(i, j) = {randomXZZMatrix(i, j), Eigen::Vector3d::Random()};
        randomYYZderivatives(i, j) = {randomYYZMatrix(i, j), Eigen::Vector3d::Random()};
        randomYZZderivatives(i, j) = {randomYZZMatrix(i, j), Eigen::Vector3d::Random()};
      }
    }
  }
};

TEST_F(AMultipoleMatrixTest, AfterResetAllMatrixEntriesAreZero) {
  dipoleMatrix.reset(arbitraryDimension);
  quadrupoleMatrix.reset(arbitraryDimension);
  octupoleMatrix.reset(arbitraryDimension);

  for (int dimension = 0; dimension < 10; ++dimension) {
    for (int i = 0; i < arbitraryDimension; ++i) {
      for (int j = 0; j < arbitraryDimension; ++j) {
        EXPECT_EQ(octupoleMatrix[dimension](i, j), 0.0);
      }
    }
  }
  for (int dimension = 0; dimension < 6; ++dimension) {
    for (int i = 0; i < arbitraryDimension; ++i) {
      for (int j = 0; j < arbitraryDimension; ++j) {
        EXPECT_EQ(quadrupoleMatrix[dimension](i, j), 0.0);
      }
    }
  }
  for (int dimension = 0; dimension < 3; ++dimension) {
    for (int i = 0; i < arbitraryDimension; ++i) {
      for (int j = 0; j < arbitraryDimension; ++j) {
        EXPECT_EQ(dipoleMatrix[dimension](i, j), 0.0);
      }
    }
  }
}

TEST_F(AMultipoleMatrixTest, MatricesAreCorrectlyAssignedWithZeroDerivative) {
  dipoleMatrix.x().get<derivOrder::zero>() = randomXMatrix;
  dipoleMatrix.y().get<derivOrder::zero>() = randomYMatrix;
  dipoleMatrix.z().get<derivOrder::zero>() = randomZMatrix;
  quadrupoleMatrix.xx().get<derivOrder::zero>() = randomXMatrix;
  quadrupoleMatrix.yy().get<derivOrder::zero>() = randomYMatrix;
  quadrupoleMatrix.zz().get<derivOrder::zero>() = randomZMatrix;
  quadrupoleMatrix.xy().get<derivOrder::zero>() = randomXYMatrix;
  quadrupoleMatrix.xz().get<derivOrder::zero>() = randomXZMatrix;
  quadrupoleMatrix.yz().get<derivOrder::zero>() = randomYZMatrix;
  octupoleMatrix.xxx().get<derivOrder::zero>() = randomXMatrix;
  octupoleMatrix.yyy().get<derivOrder::zero>() = randomYMatrix;
  octupoleMatrix.zzz().get<derivOrder::zero>() = randomZMatrix;
  octupoleMatrix.xxy().get<derivOrder::zero>() = randomXYMatrix;
  octupoleMatrix.xxz().get<derivOrder::zero>() = randomXZMatrix;
  octupoleMatrix.xyz().get<derivOrder::zero>() = randomYZMatrix;
  octupoleMatrix.xyy().get<derivOrder::zero>() = randomXYYMatrix;
  octupoleMatrix.xzz().get<derivOrder::zero>() = randomXZZMatrix;
  octupoleMatrix.yyz().get<derivOrder::zero>() = randomYYZMatrix;
  octupoleMatrix.yzz().get<derivOrder::zero>() = randomYZZMatrix;
  EXPECT_EQ(dipoleMatrix.x().get<derivOrder::zero>(), randomXMatrix);
  EXPECT_EQ(dipoleMatrix.y().get<derivOrder::zero>(), randomYMatrix);
  EXPECT_EQ(dipoleMatrix.z().get<derivOrder::zero>(), randomZMatrix);
  EXPECT_EQ(quadrupoleMatrix.xx().get<derivOrder::zero>(), randomXMatrix);
  EXPECT_EQ(quadrupoleMatrix.yy().get<derivOrder::zero>(), randomYMatrix);
  EXPECT_EQ(quadrupoleMatrix.zz().get<derivOrder::zero>(), randomZMatrix);
  EXPECT_EQ(quadrupoleMatrix.xy().get<derivOrder::zero>(), randomXYMatrix);
  EXPECT_EQ(quadrupoleMatrix.xz().get<derivOrder::zero>(), randomXZMatrix);
  EXPECT_EQ(quadrupoleMatrix.yz().get<derivOrder::zero>(), randomYZMatrix);
  EXPECT_EQ(octupoleMatrix.xxx().get<derivOrder::zero>(), randomXMatrix);
  EXPECT_EQ(octupoleMatrix.yyy().get<derivOrder::zero>(), randomYMatrix);
  EXPECT_EQ(octupoleMatrix.zzz().get<derivOrder::zero>(), randomZMatrix);
  EXPECT_EQ(octupoleMatrix.xxy().get<derivOrder::zero>(), randomXYMatrix);
  EXPECT_EQ(octupoleMatrix.xxz().get<derivOrder::zero>(), randomXZMatrix);
  EXPECT_EQ(octupoleMatrix.xyz().get<derivOrder::zero>(), randomYZMatrix);
  EXPECT_EQ(octupoleMatrix.xyy().get<derivOrder::zero>(), randomXYYMatrix);
  EXPECT_EQ(octupoleMatrix.xzz().get<derivOrder::zero>(), randomXZZMatrix);
  EXPECT_EQ(octupoleMatrix.yyz().get<derivOrder::zero>(), randomYYZMatrix);
  EXPECT_EQ(octupoleMatrix.yzz().get<derivOrder::zero>(), randomYZZMatrix);
}

TEST_F(AMultipoleMatrixTest, MatricesAreCorrectlyAssignedWithZeroDerivativeAndBlockMethod) {
  dipoleMatrix.reset(randomXMatrix.rows());
  dipoleMatrix.x().get<derivOrder::zero>().block(0, 0, 5, 5) = randomXMatrix.block(0, 0, 5, 5);

  EXPECT_EQ(dipoleMatrix.x().get<derivOrder::zero>().block(0, 0, 5, 5), randomXMatrix.block(0, 0, 5, 5));
}

TEST_F(AMultipoleMatrixTest, MatricesAreCorrectlyAssignedWithDerivatives) {
  dipoleMatrix.x().get<derivOrder::one>() = randomXderivatives;
  dipoleMatrix.y().get<derivOrder::one>() = randomYderivatives;
  dipoleMatrix.z().get<derivOrder::one>() = randomZderivatives;
  quadrupoleMatrix.xx().get<derivOrder::one>() = randomXderivatives;
  quadrupoleMatrix.yy().get<derivOrder::one>() = randomYderivatives;
  quadrupoleMatrix.zz().get<derivOrder::one>() = randomZderivatives;
  quadrupoleMatrix.xy().get<derivOrder::one>() = randomXYderivatives;
  quadrupoleMatrix.xz().get<derivOrder::one>() = randomXZderivatives;
  quadrupoleMatrix.yz().get<derivOrder::one>() = randomYZderivatives;
  octupoleMatrix.xxx().get<derivOrder::one>() = randomXderivatives;
  octupoleMatrix.yyy().get<derivOrder::one>() = randomYderivatives;
  octupoleMatrix.zzz().get<derivOrder::one>() = randomZderivatives;
  octupoleMatrix.xxy().get<derivOrder::one>() = randomXYderivatives;
  octupoleMatrix.xxz().get<derivOrder::one>() = randomXZderivatives;
  octupoleMatrix.xyz().get<derivOrder::one>() = randomYZderivatives;
  octupoleMatrix.xyy().get<derivOrder::one>() = randomXYYderivatives;
  octupoleMatrix.xzz().get<derivOrder::one>() = randomXZZderivatives;
  octupoleMatrix.yyz().get<derivOrder::one>() = randomYYZderivatives;
  octupoleMatrix.yzz().get<derivOrder::one>() = randomYZZderivatives;
  for (int i = 0; i < arbitraryDimension; ++i) {
    for (int j = 0; j < arbitraryDimension; ++j) {
      EXPECT_EQ(dipoleMatrix.x().get<derivOrder::one>()(i, j).derivatives(), randomXderivatives(i, j).derivatives());
      EXPECT_EQ(dipoleMatrix.y().get<derivOrder::one>()(i, j).derivatives(), randomYderivatives(i, j).derivatives());
      EXPECT_EQ(dipoleMatrix.z().get<derivOrder::one>()(i, j).derivatives(), randomZderivatives(i, j).derivatives());
      EXPECT_EQ(quadrupoleMatrix.xx().get<derivOrder::one>()(i, j).derivatives(), randomXderivatives(i, j).derivatives());
      EXPECT_EQ(quadrupoleMatrix.yy().get<derivOrder::one>()(i, j).derivatives(), randomYderivatives(i, j).derivatives());
      EXPECT_EQ(quadrupoleMatrix.zz().get<derivOrder::one>()(i, j).derivatives(), randomZderivatives(i, j).derivatives());
      EXPECT_EQ(quadrupoleMatrix.xy().get<derivOrder::one>()(i, j).derivatives(), randomXYderivatives(i, j).derivatives());
      EXPECT_EQ(quadrupoleMatrix.xz().get<derivOrder::one>()(i, j).derivatives(), randomXZderivatives(i, j).derivatives());
      EXPECT_EQ(quadrupoleMatrix.yz().get<derivOrder::one>()(i, j).derivatives(), randomYZderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.xxx().get<derivOrder::one>()(i, j).derivatives(), randomXderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.yyy().get<derivOrder::one>()(i, j).derivatives(), randomYderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.zzz().get<derivOrder::one>()(i, j).derivatives(), randomZderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.xxy().get<derivOrder::one>()(i, j).derivatives(), randomXYderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.xxz().get<derivOrder::one>()(i, j).derivatives(), randomXZderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.xyz().get<derivOrder::one>()(i, j).derivatives(), randomYZderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.xyy().get<derivOrder::one>()(i, j).derivatives(), randomXYYderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.xzz().get<derivOrder::one>()(i, j).derivatives(), randomXZZderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.yyz().get<derivOrder::one>()(i, j).derivatives(), randomYYZderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.yzz().get<derivOrder::one>()(i, j).derivatives(), randomYZZderivatives(i, j).derivatives());
    }
  }
}

TEST_F(AMultipoleMatrixTest, ThrowsOnlyIfDimensionOutOfBounds) {
  EXPECT_THROW(dipoleMatrix[4], std::exception);
  EXPECT_THROW(dipoleMatrix[-2], std::exception);
  EXPECT_NO_THROW(dipoleMatrix[0]);
  EXPECT_NO_THROW(dipoleMatrix[1]);
  EXPECT_NO_THROW(dipoleMatrix[2]);
  EXPECT_THROW(quadrupoleMatrix[6], std::exception);
  EXPECT_THROW(quadrupoleMatrix[-2], std::exception);
  EXPECT_NO_THROW(quadrupoleMatrix[0]);
  EXPECT_NO_THROW(quadrupoleMatrix[1]);
  EXPECT_NO_THROW(quadrupoleMatrix[2]);
  EXPECT_NO_THROW(quadrupoleMatrix[3]);
  EXPECT_NO_THROW(quadrupoleMatrix[4]);
  EXPECT_NO_THROW(quadrupoleMatrix[5]);
}
} // namespace Tests
} // namespace Utils
} // namespace Scine
