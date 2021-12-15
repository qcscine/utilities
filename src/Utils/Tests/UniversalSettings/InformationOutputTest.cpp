/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/UniversalSettings/DescriptorCollection.h>
#include <Utils/UniversalSettings/InformationOutput.h>
#include <Utils/UniversalSettings/SettingPopulator.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <gmock/gmock.h>
#include <boost/dll.hpp>
#include <boost/filesystem/path.hpp>
#include <fstream>
#include <regex>
#include <string>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class InformationOutputTest InformationOutputTest.cpp
 * @brief Comprises tests for the class Scine::Utils::UniversalSettings::InformationOutput.
 * @test
 */
class InformationOutputTest : public Test {
 public:
  boost::filesystem::path dummyFile;

  std::ofstream openDummyFile() {
    dummyFile = boost::dll::program_location().parent_path();
    dummyFile /= "test_print";
    std::ofstream fout;
    fout.open(dummyFile.string());
    if (!fout.is_open()) {
      throw std::runtime_error("Problem when opening/creating file " + dummyFile.string());
    }
    return fout;
  }

  std::string readDummyFile() {
    std::ifstream fin;
    fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    fin.open(dummyFile.string());
    auto content = std::string(std::istreambuf_iterator<char>{fin}, {});
    fin.close();
    int i = std::remove(dummyFile.c_str());
    if (i != 0) {
      throw std::runtime_error("Could not delete file " + dummyFile.string());
    }
    return content;
  }
};

TEST_F(InformationOutputTest, CanPrintTitle) {
  using namespace UniversalSettings;
  auto fields = DescriptorCollection("TestDescription");
  ASSERT_TRUE(fields.empty());
  std::string key = "key";
  auto fout = openDummyFile();
  InformationOutput::print(key, fields, fout, 1, true);
  fout.close();
  auto content = readDummyFile();
  std::regex regex("TestDescription");
  std::smatch smatch;
  ASSERT_TRUE(std::regex_search(content, smatch, regex));
  fout = openDummyFile();
  InformationOutput::printLong(key, fields, fout, 2);
  fout.close();
  content = readDummyFile();
  ASSERT_TRUE(std::regex_search(content, smatch, regex));
}

TEST_F(InformationOutputTest, CanPrintSettings) {
  using namespace UniversalSettings;
  auto fields = DescriptorCollection("TestDescription");
  ASSERT_TRUE(fields.empty());
  SettingPopulator::populateLcaoSettings(fields);
  std::string key = "key";
  auto fout = openDummyFile();
  InformationOutput::print(key, fields, fout, 1, true);
  fout.close();
  auto content = readDummyFile();
  std::vector<std::string> settings = {SettingsNames::spinMode, SettingsNames::spinMultiplicity,
                                       SettingsNames::molecularCharge, SettingsNames::temperature};
  for (const auto& setting : settings) {
    std::regex regex(setting);
    std::smatch smatch;
    ASSERT_TRUE(std::regex_search(content, smatch, regex));
  }
  fout = openDummyFile();
  InformationOutput::printLong(key, fields, fout, 2);
  fout.close();
  content = readDummyFile();
  for (const auto& setting : settings) {
    std::regex regex(setting);
    std::smatch smatch;
    ASSERT_TRUE(std::regex_search(content, smatch, regex));
  }
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
