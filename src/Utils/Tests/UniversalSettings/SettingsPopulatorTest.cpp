/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/UniversalSettings/DescriptorCollection.h>
#include <Utils/UniversalSettings/SettingPopulator.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <gmock/gmock.h>
#include <string>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class SettingsPopulatorTest SettingsPopulatorTest.cpp
 * @brief Comprises tests for the class Scine::Utils::UniversalSettings::SettingsPopulator.
 * @test
 */
class SettingsPopulatorTest : public Test {};

TEST_F(SettingsPopulatorTest, CanAddToFields) {
  using namespace UniversalSettings;
  auto fields = DescriptorCollection();
  ASSERT_TRUE(fields.empty());
  SettingPopulator::addLogOption(fields);
  SettingPopulator::populateLcaoSettings(fields);
  SettingPopulator::populateScfSettings(fields);
  SettingPopulator::populateSemiEmpiricalSettings(fields, "test.param");
  ASSERT_TRUE(fields.exists(SettingsNames::loggerVerbosity));
  ASSERT_TRUE(fields.exists(SettingsNames::molecularCharge));
  ASSERT_TRUE(fields.exists(SettingsNames::spinMultiplicity));
  ASSERT_TRUE(fields.exists(SettingsNames::spinMode));
  ASSERT_TRUE(fields.exists(SettingsNames::selfConsistenceCriterion));
  ASSERT_TRUE(fields.exists(SettingsNames::maxScfIterations));
  ASSERT_TRUE(fields.exists(SettingsNames::mixer));
  ASSERT_TRUE(fields.exists(SettingsNames::temperature));
  ASSERT_TRUE(fields.exists(SettingsNames::electronicTemperature));
  ASSERT_TRUE(fields.exists(SettingsNames::symmetryNumber));
  ASSERT_TRUE(fields.exists(SettingsNames::methodParameters));
  GenericValue param = fields.get(SettingsNames::methodParameters).getDefaultValue();
  std::string p = param;
  ASSERT_EQ(p, "test.param");
}

TEST_F(SettingsPopulatorTest, ScfMixerIsSelfConsistent) {
  using namespace UniversalSettings;
  ASSERT_EQ(SettingPopulator::stringToScfMixer(SettingPopulator::scfMixerToString(scf_mixer_t::fock_diis)),
            scf_mixer_t::fock_diis);
  ASSERT_EQ(SettingPopulator::stringToScfMixer(SettingPopulator::scfMixerToString(scf_mixer_t::ediis)), scf_mixer_t::ediis);
  ASSERT_EQ(SettingPopulator::stringToScfMixer(SettingPopulator::scfMixerToString(scf_mixer_t::ediis_diis)),
            scf_mixer_t::ediis_diis);
  ASSERT_EQ(SettingPopulator::stringToScfMixer(SettingPopulator::scfMixerToString(scf_mixer_t::none)), scf_mixer_t::none);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
