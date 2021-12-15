/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/TimeDependent/TransitionDipoleCalculator.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {
TEST(TransitionDipoleMomentTest, CanTransformTransitionDipoleToOscillatorStrength) {
  Eigen::Matrix3Xd transDipMat(3, 2);
  transDipMat << 1.0, -1.0, 0.0, 0.0, 0.0, 0.0;
  Eigen::VectorXd eigenvalues = Eigen::Vector2d(1., 2.);

  Eigen::VectorXd oscStrengths =
      TransitionDipoleCalculator::transitionDipoleMomentToOscillatorStrength(transDipMat, eigenvalues);

  ASSERT_EQ(oscStrengths.size(), 2);
  Eigen::VectorXd actualOscStrengths(2);
  actualOscStrengths << 2.0 / 3.0, 4.0 / 3.0;

  ASSERT_DOUBLE_EQ(oscStrengths(0), actualOscStrengths(0));
  ASSERT_DOUBLE_EQ(oscStrengths(1), actualOscStrengths(1));
}
} // namespace Tests
} // namespace Utils
} // namespace Scine
