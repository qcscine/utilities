/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
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
  dipoleMatrix.x().get<DerivativeOrder::Zero>() = randomXMatrix;
  dipoleMatrix.y().get<DerivativeOrder::Zero>() = randomYMatrix;
  dipoleMatrix.z().get<DerivativeOrder::Zero>() = randomZMatrix;
  quadrupoleMatrix.xx().get<DerivativeOrder::Zero>() = randomXMatrix;
  quadrupoleMatrix.yy().get<DerivativeOrder::Zero>() = randomYMatrix;
  quadrupoleMatrix.zz().get<DerivativeOrder::Zero>() = randomZMatrix;
  quadrupoleMatrix.xy().get<DerivativeOrder::Zero>() = randomXYMatrix;
  quadrupoleMatrix.xz().get<DerivativeOrder::Zero>() = randomXZMatrix;
  quadrupoleMatrix.yz().get<DerivativeOrder::Zero>() = randomYZMatrix;
  octupoleMatrix.xxx().get<DerivativeOrder::Zero>() = randomXMatrix;
  octupoleMatrix.yyy().get<DerivativeOrder::Zero>() = randomYMatrix;
  octupoleMatrix.zzz().get<DerivativeOrder::Zero>() = randomZMatrix;
  octupoleMatrix.xxy().get<DerivativeOrder::Zero>() = randomXYMatrix;
  octupoleMatrix.xxz().get<DerivativeOrder::Zero>() = randomXZMatrix;
  octupoleMatrix.xyz().get<DerivativeOrder::Zero>() = randomYZMatrix;
  octupoleMatrix.xyy().get<DerivativeOrder::Zero>() = randomXYYMatrix;
  octupoleMatrix.xzz().get<DerivativeOrder::Zero>() = randomXZZMatrix;
  octupoleMatrix.yyz().get<DerivativeOrder::Zero>() = randomYYZMatrix;
  octupoleMatrix.yzz().get<DerivativeOrder::Zero>() = randomYZZMatrix;
  EXPECT_EQ(dipoleMatrix.x().get<DerivativeOrder::Zero>(), randomXMatrix);
  EXPECT_EQ(dipoleMatrix.y().get<DerivativeOrder::Zero>(), randomYMatrix);
  EXPECT_EQ(dipoleMatrix.z().get<DerivativeOrder::Zero>(), randomZMatrix);
  EXPECT_EQ(quadrupoleMatrix.xx().get<DerivativeOrder::Zero>(), randomXMatrix);
  EXPECT_EQ(quadrupoleMatrix.yy().get<DerivativeOrder::Zero>(), randomYMatrix);
  EXPECT_EQ(quadrupoleMatrix.zz().get<DerivativeOrder::Zero>(), randomZMatrix);
  EXPECT_EQ(quadrupoleMatrix.xy().get<DerivativeOrder::Zero>(), randomXYMatrix);
  EXPECT_EQ(quadrupoleMatrix.xz().get<DerivativeOrder::Zero>(), randomXZMatrix);
  EXPECT_EQ(quadrupoleMatrix.yz().get<DerivativeOrder::Zero>(), randomYZMatrix);
  EXPECT_EQ(octupoleMatrix.xxx().get<DerivativeOrder::Zero>(), randomXMatrix);
  EXPECT_EQ(octupoleMatrix.yyy().get<DerivativeOrder::Zero>(), randomYMatrix);
  EXPECT_EQ(octupoleMatrix.zzz().get<DerivativeOrder::Zero>(), randomZMatrix);
  EXPECT_EQ(octupoleMatrix.xxy().get<DerivativeOrder::Zero>(), randomXYMatrix);
  EXPECT_EQ(octupoleMatrix.xxz().get<DerivativeOrder::Zero>(), randomXZMatrix);
  EXPECT_EQ(octupoleMatrix.xyz().get<DerivativeOrder::Zero>(), randomYZMatrix);
  EXPECT_EQ(octupoleMatrix.xyy().get<DerivativeOrder::Zero>(), randomXYYMatrix);
  EXPECT_EQ(octupoleMatrix.xzz().get<DerivativeOrder::Zero>(), randomXZZMatrix);
  EXPECT_EQ(octupoleMatrix.yyz().get<DerivativeOrder::Zero>(), randomYYZMatrix);
  EXPECT_EQ(octupoleMatrix.yzz().get<DerivativeOrder::Zero>(), randomYZZMatrix);
}

TEST_F(AMultipoleMatrixTest, MatricesAreCorrectlyAssignedWithZeroDerivativeAndBlockMethod) {
  dipoleMatrix.reset(randomXMatrix.rows());
  dipoleMatrix.x().get<DerivativeOrder::Zero>().block(0, 0, 5, 5) = randomXMatrix.block(0, 0, 5, 5);

  EXPECT_EQ(dipoleMatrix.x().get<DerivativeOrder::Zero>().block(0, 0, 5, 5), randomXMatrix.block(0, 0, 5, 5));
}

TEST_F(AMultipoleMatrixTest, MatricesAreCorrectlyAssignedWithDerivatives) {
  dipoleMatrix.x().get<DerivativeOrder::One>() = randomXderivatives;
  dipoleMatrix.y().get<DerivativeOrder::One>() = randomYderivatives;
  dipoleMatrix.z().get<DerivativeOrder::One>() = randomZderivatives;
  quadrupoleMatrix.xx().get<DerivativeOrder::One>() = randomXderivatives;
  quadrupoleMatrix.yy().get<DerivativeOrder::One>() = randomYderivatives;
  quadrupoleMatrix.zz().get<DerivativeOrder::One>() = randomZderivatives;
  quadrupoleMatrix.xy().get<DerivativeOrder::One>() = randomXYderivatives;
  quadrupoleMatrix.xz().get<DerivativeOrder::One>() = randomXZderivatives;
  quadrupoleMatrix.yz().get<DerivativeOrder::One>() = randomYZderivatives;
  octupoleMatrix.xxx().get<DerivativeOrder::One>() = randomXderivatives;
  octupoleMatrix.yyy().get<DerivativeOrder::One>() = randomYderivatives;
  octupoleMatrix.zzz().get<DerivativeOrder::One>() = randomZderivatives;
  octupoleMatrix.xxy().get<DerivativeOrder::One>() = randomXYderivatives;
  octupoleMatrix.xxz().get<DerivativeOrder::One>() = randomXZderivatives;
  octupoleMatrix.xyz().get<DerivativeOrder::One>() = randomYZderivatives;
  octupoleMatrix.xyy().get<DerivativeOrder::One>() = randomXYYderivatives;
  octupoleMatrix.xzz().get<DerivativeOrder::One>() = randomXZZderivatives;
  octupoleMatrix.yyz().get<DerivativeOrder::One>() = randomYYZderivatives;
  octupoleMatrix.yzz().get<DerivativeOrder::One>() = randomYZZderivatives;
  for (int i = 0; i < arbitraryDimension; ++i) {
    for (int j = 0; j < arbitraryDimension; ++j) {
      EXPECT_EQ(dipoleMatrix.x().get<DerivativeOrder::One>()(i, j).derivatives(), randomXderivatives(i, j).derivatives());
      EXPECT_EQ(dipoleMatrix.y().get<DerivativeOrder::One>()(i, j).derivatives(), randomYderivatives(i, j).derivatives());
      EXPECT_EQ(dipoleMatrix.z().get<DerivativeOrder::One>()(i, j).derivatives(), randomZderivatives(i, j).derivatives());
      EXPECT_EQ(quadrupoleMatrix.xx().get<DerivativeOrder::One>()(i, j).derivatives(), randomXderivatives(i, j).derivatives());
      EXPECT_EQ(quadrupoleMatrix.yy().get<DerivativeOrder::One>()(i, j).derivatives(), randomYderivatives(i, j).derivatives());
      EXPECT_EQ(quadrupoleMatrix.zz().get<DerivativeOrder::One>()(i, j).derivatives(), randomZderivatives(i, j).derivatives());
      EXPECT_EQ(quadrupoleMatrix.xy().get<DerivativeOrder::One>()(i, j).derivatives(),
                randomXYderivatives(i, j).derivatives());
      EXPECT_EQ(quadrupoleMatrix.xz().get<DerivativeOrder::One>()(i, j).derivatives(),
                randomXZderivatives(i, j).derivatives());
      EXPECT_EQ(quadrupoleMatrix.yz().get<DerivativeOrder::One>()(i, j).derivatives(),
                randomYZderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.xxx().get<DerivativeOrder::One>()(i, j).derivatives(), randomXderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.yyy().get<DerivativeOrder::One>()(i, j).derivatives(), randomYderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.zzz().get<DerivativeOrder::One>()(i, j).derivatives(), randomZderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.xxy().get<DerivativeOrder::One>()(i, j).derivatives(), randomXYderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.xxz().get<DerivativeOrder::One>()(i, j).derivatives(), randomXZderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.xyz().get<DerivativeOrder::One>()(i, j).derivatives(), randomYZderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.xyy().get<DerivativeOrder::One>()(i, j).derivatives(),
                randomXYYderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.xzz().get<DerivativeOrder::One>()(i, j).derivatives(),
                randomXZZderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.yyz().get<DerivativeOrder::One>()(i, j).derivatives(),
                randomYYZderivatives(i, j).derivatives());
      EXPECT_EQ(octupoleMatrix.yzz().get<DerivativeOrder::One>()(i, j).derivatives(),
                randomYZZderivatives(i, j).derivatives());
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
