/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Math/AtomicSecondDerivativeCollection.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class ASecondDerivativeCollectionTest SecondDerivativeCollectionTest.cpp
 * @brief Comprises tests for the class Scine::Utils::SecondDerivativeCollection.
 * @test
 */
class ASecondDerivativeCollectionTest : public Test {};

TEST_F(ASecondDerivativeCollectionTest, WorksCorrectly) {
  AtomicSecondDerivativeCollection asdc1 = AtomicSecondDerivativeCollection();
  ASSERT_EQ(asdc1.size(), 0);
  ASSERT_TRUE(asdc1.empty());

  AtomicSecondDerivativeCollection asdc2 = AtomicSecondDerivativeCollection(2);
  ASSERT_EQ(asdc2.size(), 2);
  asdc2.resize(1);
  ASSERT_EQ(asdc2.size(), 1);

  auto ahs = asdc2.getAtomicHessians();
  ASSERT_EQ(ahs.size(), 1);

  asdc2.setZero();
  Eigen::Matrix3d ah = asdc2.getAtomicHessian(0);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_NEAR(ah(i, j), 0.0, 1.0e-10);
    }
  }

  asdc2.clear();
  ASSERT_EQ(asdc2.size(), 0);
}

} // namespace Tests
} // namespace Utils
} // namespace Scine