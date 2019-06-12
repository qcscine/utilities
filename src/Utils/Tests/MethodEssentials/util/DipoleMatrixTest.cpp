/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/MethodEssentials/util/DipoleMatrix.h>
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Utils {
namespace Tests {

class ADipoleMatrixTest : public Test {
 public:
  int arbitraryDimension = 10;
  DipoleMatrix dipoleMatrix;
  Eigen::MatrixXd randomXMatrix, randomYMatrix, randomZMatrix;

 private:
  void SetUp() final {
    randomXMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    randomYMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    randomZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
  }
};

TEST_F(ADipoleMatrixTest, AfterResetAllMatrixEntriesAreZero) {
  dipoleMatrix.reset(arbitraryDimension);

  for (int dimension = 0; dimension < 3; ++dimension) {
    for (int i = 0; i < arbitraryDimension; ++i) {
      for (int j = 0; j < arbitraryDimension; ++j) {
        ASSERT_EQ(dipoleMatrix[dimension](i, j), 0.0);
      }
    }
  }
}

TEST_F(ADipoleMatrixTest, MatricesAreCorrectlyAssigned) {
  dipoleMatrix.x() = randomXMatrix;
  dipoleMatrix.y() = randomYMatrix;
  dipoleMatrix.z() = randomZMatrix;

  for (int i = 0; i < arbitraryDimension; ++i) {
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_EQ(dipoleMatrix.x()(i, j), randomXMatrix(i, j));
    }
  }
  for (int i = 0; i < arbitraryDimension; ++i) {
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_EQ(dipoleMatrix.y()(i, j), randomYMatrix(i, j));
    }
  }
  for (int i = 0; i < arbitraryDimension; ++i) {
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_EQ(dipoleMatrix.z()(i, j), randomZMatrix(i, j));
    }
  }
}

TEST_F(ADipoleMatrixTest, ThrowsOnlyIfDimensionOutOfBounds) {
  ASSERT_THROW(dipoleMatrix[4], std::exception);
  ASSERT_THROW(dipoleMatrix[-2], std::exception);
  ASSERT_NO_THROW(dipoleMatrix[0]);
  ASSERT_NO_THROW(dipoleMatrix[1]);
  ASSERT_NO_THROW(dipoleMatrix[2]);
}
} // namespace Tests
} // namespace Utils
} // namespace Scine
