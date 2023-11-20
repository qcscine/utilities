/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/UniversalSettings/BoolDescriptor.h>
#include <gmock/gmock.h>
#include <string>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class BoolDescriptorTest BoolDescriptorTest.cpp
 * @brief Comprises tests for the class Scine::Utils::UniversalSettings::BoolDescriptor.
 * @test
 */
class BoolDescriptorTest : public Test {};

TEST_F(BoolDescriptorTest, ReturnsCorrectPropertyDescription) {
  UniversalSettings::BoolDescriptor boolDescriptor("Initial description");
  ASSERT_EQ(boolDescriptor.getPropertyDescription(), "Initial description");
  boolDescriptor.setPropertyDescription("Second description");
  ASSERT_EQ(boolDescriptor.getPropertyDescription(), "Second description");
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
