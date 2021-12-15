/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/LennardJonesCalculator/LennardJonesCalculator.h>
#include <Utils/LennardJonesCalculator/LennardJonesCalculatorSettings.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

class ALennardJonesTest : public Test {
 public:
  LennardJonesCalculator calculator;

 private:
  void SetUp() final {
  }
};

TEST_F(ALennardJonesTest, SettingsAreSetCorrectly) {
  calculator.settings().modifyDouble(Utils::SettingsNames::lennardJonesSigma, 4.0);
  calculator.settings().modifyDouble(Utils::SettingsNames::lennardJonesEpsilon, 3.0);
  calculator.settings().modifyDouble(Utils::SettingsNames::lennardJonesCutoff, 2.0);
  calculator.settings().modifyBool(Utils::SettingsNames::lennardJonesUsePBCs, true);
  calculator.settings().modifyDouble(Utils::SettingsNames::lennardJonesBoxsize, 5.0);

  ASSERT_THAT(calculator.settings().getDouble(Utils::SettingsNames::lennardJonesSigma), Eq(4.0));
  ASSERT_THAT(calculator.settings().getDouble(Utils::SettingsNames::lennardJonesEpsilon), Eq(3.0));
  ASSERT_THAT(calculator.settings().getDouble(Utils::SettingsNames::lennardJonesCutoff), Eq(2.0));
  ASSERT_TRUE(calculator.settings().getBool(Utils::SettingsNames::lennardJonesUsePBCs));
  ASSERT_THAT(calculator.settings().getDouble(Utils::SettingsNames::lennardJonesBoxsize), Eq(5.0));
}

TEST_F(ALennardJonesTest, SettingsCompatibilityIsChecked) {
  calculator.settings().modifyDouble(Utils::SettingsNames::lennardJonesCutoff, 4.0);
  calculator.settings().modifyBool(Utils::SettingsNames::lennardJonesUsePBCs, true);
  calculator.settings().modifyDouble(Utils::SettingsNames::lennardJonesBoxsize, 5.0);

  std::stringstream stream("2\n\n"
                           "Ar    0.00000000   0.00000000   0.00000000\n"
                           "Ar    0.00000000   0.00000000   1.50000000\n");

  auto structure = Utils::XyzStreamHandler::read(stream);

  // Should throw due to incompatible box size and cut off
  ASSERT_THROW(calculator.setStructure(structure), Core::InitializationException);
}

TEST_F(ALennardJonesTest, CanCalculateWithoutPBCs) {
  std::stringstream stream("3\n\n"
                           "Ar    0.80000000   0.00000000   0.00000000\n"
                           "Ar    0.00000000   0.34000000   2.00000000\n"
                           "Ar    0.00000000   0.00000000   4.50000000\n");

  auto structure = Utils::XyzStreamHandler::read(stream);
  calculator.setStructure(structure);
  calculator.calculate("test description");
  ASSERT_THAT(calculator.results().get<Property::Description>(), Eq("test description"));
  // Reference energy and gradients obtained from CP2K version 7.0, git:a515b01
  EXPECT_NEAR(calculator.results().get<Property::Energy>(), 0.320842529947627, 1e-4);
  Eigen::Matrix3d cp2kForces;
  cp2kForces << 0.30823624, -0.13100910, -0.77052665, -0.30825671, 0.14714612, 0.65198721, 0.00002047, -0.01613702, 0.11853944;
  ASSERT_TRUE(calculator.results().get<Property::Gradients>().isApprox(-1.0 * cp2kForces, 1e-4));

  // With a cut-off smaller than the interatomic distances the energy should be zero
  calculator.settings().modifyDouble(Utils::SettingsNames::lennardJonesCutoff, 0.2);
  calculator.calculate("");
  EXPECT_NEAR(calculator.results().get<Property::Energy>(), 0.0, 1e-4);
  ASSERT_TRUE(calculator.results().get<Property::Gradients>().isZero());
}

TEST_F(ALennardJonesTest, CanCalculateWithPBCs) {
  std::stringstream stream("4\n\n"
                           "Ar    0.80000000   0.00000000  42.00000000\n"
                           "Ar    0.00000000   0.34000000  19.00000000\n"
                           "Ar    0.10000000   0.00000000  23.00000000\n"
                           "Ar   21.00000000  11.00000000   0.10000000\n");

  auto structure = Utils::XyzStreamHandler::read(stream);
  calculator.setStructure(structure);
  calculator.settings().modifyBool(Utils::SettingsNames::lennardJonesUsePBCs, true);
  calculator.calculate("test description");
  ASSERT_THAT(calculator.results().get<Property::Description>(), Eq("test description"));
  // Reference energy and gradients obtained from CP2K version 7.0, git:a515b01
  EXPECT_NEAR(calculator.results().get<Property::Energy>(), 0.840893103151347, 1e-4);
  Eigen::Matrix<double, 4, 3> cp2kForces;
  cp2kForces << 1.00491747, -0.36529297, 1.51911449, -0.85950985, 0.36528270, -1.96914319, -0.14540762, 0.00001026,
      0.45002870, 0.00000000, 0.00000000, 0.00000000;
  ASSERT_TRUE(calculator.results().get<Property::Gradients>().isApprox(-1.0 * cp2kForces, 1e-4));

  // With a cut-off smaller than the interatomic distances the energy should be zero
  calculator.settings().modifyDouble(Utils::SettingsNames::lennardJonesCutoff, 0.2);
  calculator.calculate("");
  EXPECT_NEAR(calculator.results().get<Property::Energy>(), 0.0, 1e-4);
  ASSERT_TRUE(calculator.results().get<Property::Gradients>().isZero());
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
