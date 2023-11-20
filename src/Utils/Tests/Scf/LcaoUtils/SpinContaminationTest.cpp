/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Scf/LcaoUtils/SpinContamination.h"
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace LcaoUtils {
namespace Tests {

/**
 * @class Scine::Utils::Tests::SpinContamination
 * @brief Very small test for the evaluation of the S^2 value.
 * @test
 */
TEST(SpinContaminationTests, evaluteS2correctly) {
  Eigen::Matrix2d coeffAlpha;
  coeffAlpha << 0.5445586948, 1.2620659394, 0.5445586948, -1.2620659394;
  Eigen::Matrix2d coeffBeta = coeffAlpha;
  Eigen::Matrix2d overlap;
  overlap << 1.0000000000, 0.6860894499, 0.6860894499, 1.0000000000;

  auto s2 = S2(overlap, coeffAlpha, coeffBeta, 1, 1);

  EXPECT_THAT(s2.first, DoubleNear(0.0, 1e-10));
  EXPECT_THAT(s2.second, DoubleNear(0.0, 1e-10));
}

} /* namespace Tests */
} /* namespace LcaoUtils */
} /* namespace Utils */
} /* namespace Scine */
