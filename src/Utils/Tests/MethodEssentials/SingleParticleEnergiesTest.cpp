/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/MethodEssentials/util/SingleParticleEnergies.h>
#include <gmock/gmock.h>
#include <Eigen/Core>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

class ASingleParticleEnergies : public Test {
 public:
  SingleParticleEnergies e;
};

TEST_F(ASingleParticleEnergies, InitializesToZero) {
  ASSERT_THAT(e.getRestrictedNLevels(), Eq(0));
}

TEST_F(ASingleParticleEnergies, HasCorrectNumberOfLevels) {
  Eigen::VectorXd values(4);
  values << -2.2, 0, 3.3, 9.9;

  e.setRestricted(values);

  ASSERT_THAT(e.getRestrictedNLevels(), Eq(4));
}

TEST_F(ASingleParticleEnergies, ReturnsLevelEnergy) {
  Eigen::VectorXd values(4);
  values << -2.2, 0, 3.3, 9.9;

  e.setRestricted(values);

  ASSERT_THAT(e.getRestrictedLevelEnergy(2), DoubleEq(3.3));
}

TEST_F(ASingleParticleEnergies, ReturnsEnergiesVector) {
  Eigen::VectorXd values(4);
  values << -2.2, 0, 3.3, 9.9;

  e.setRestricted(values);

  for (int i = 0; i < 4; i++) {
    ASSERT_THAT(e.getRestrictedLevelEnergy(i), DoubleEq(values(i)));
  }
}

TEST_F(ASingleParticleEnergies, CanSetUnrestrictedEnergies) {
  Eigen::VectorXd valuesAlpha(5), valuesBeta(5);
  valuesAlpha << -2.2, 0, 3.3, 9.9, 10;
  valuesBeta << -3.2, -1, 1.3, 2.9, 20;

  e.setUnrestricted(valuesAlpha, valuesBeta);

  ASSERT_THAT(e.getAlphaLevelEnergy(4), DoubleEq(10));
  ASSERT_THAT(e.getBetaLevelEnergy(3), DoubleEq(2.9));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
