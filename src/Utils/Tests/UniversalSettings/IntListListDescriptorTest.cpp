/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/UniversalSettings/IntListListDescriptor.h>
#include <gmock/gmock.h>
#include <string>
#include <vector>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class BoolDescriptorTest BoolDescriptorTest.cpp
 * @brief Comprises tests for the class Scine::Utils::UniversalSettings::BoolDescriptor.
 * @test
 */
class IntListListDescriptorTest : public Test {};

TEST_F(IntListListDescriptorTest, valueHandlingAndDescription) {
  UniversalSettings::IntListListDescriptor intListListDescriptor("Initial description");
  ASSERT_EQ(intListListDescriptor.getPropertyDescription(), "Initial description");
  UniversalSettings::IntListListDescriptor::IntListList defaultValue = {{0, 1, 2}, {3, 4, 5}};
  intListListDescriptor.setDefaultValue(defaultValue);
  auto value = intListListDescriptor.getDefaultValue();
  for (unsigned int i = 0; i < value.size(); ++i) {
    for (unsigned int j = 0; j < value[i].size(); ++j) {
      EXPECT_EQ(value[i][j], defaultValue[i][j]);
    }
  }
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */