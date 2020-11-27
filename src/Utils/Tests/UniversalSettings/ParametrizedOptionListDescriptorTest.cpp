/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/UniversalSettings/DescriptorCollection.h"
#include <Utils/UniversalSettings/ParametrizedOptionListDescriptor.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class ParametrizedOptionListDescriptorTest ParametrizedOptionListDescriptorTest.cpp
 * @brief Comprises tests for the class Scine::Utils::UniversalSettings::ParametrizedOptionListDescriptor.
 * @test
 */
class ParametrizedOptionListDescriptorTest : public Test {};

TEST_F(ParametrizedOptionListDescriptorTest, ReturnsCorrectDescription) {
  auto descriptor = UniversalSettings::ParametrizedOptionListDescriptor("Dummy property description");
  ASSERT_EQ(descriptor.getPropertyDescription(), "Dummy property description");
  descriptor.setPropertyDescription("New property description");
  ASSERT_EQ(descriptor.getPropertyDescription(), "New property description");
}

TEST_F(ParametrizedOptionListDescriptorTest, ReturnsCorrectOptionCountAfterInitialization) {
  auto descriptor = UniversalSettings::ParametrizedOptionListDescriptor("Dummy property description");
  ASSERT_EQ(descriptor.optionCount(), 0);
}

TEST_F(ParametrizedOptionListDescriptorTest, CanAddOptions) {
  auto descriptor = UniversalSettings::ParametrizedOptionListDescriptor("Dummy property description");
  ASSERT_FALSE(descriptor.optionExists("New dummy option"));
  descriptor.addOption("New dummy option");
  ASSERT_EQ(descriptor.optionCount(), 1);
  ASSERT_TRUE(descriptor.optionExists("New dummy option"));

  ASSERT_FALSE(descriptor.optionExists("Second new option"));
  auto optionSettings = UniversalSettings::DescriptorCollection("New descriptor collection");
  descriptor.addOption("Second new option", optionSettings);
  ASSERT_EQ(descriptor.optionCount(), 2);
  ASSERT_TRUE(descriptor.optionExists("Second new option"));
}

TEST_F(ParametrizedOptionListDescriptorTest, CanSetDefaultOption) {
  auto descriptor = UniversalSettings::ParametrizedOptionListDescriptor("Dummy property description");
  ASSERT_THROW(descriptor.getDefaultIndex(), UniversalSettings::EmptyOptionListException);
  ASSERT_THROW(descriptor.getDefaultOption(), UniversalSettings::EmptyOptionListException);
  ASSERT_THROW(descriptor.setDefaultOption("Nonexisting option"), UniversalSettings::OptionDoesNotExistException);

  descriptor.addOption("New dummy option");
  descriptor.addOption("Second new option");
  ASSERT_EQ(descriptor.getDefaultIndex(), 0);
  ASSERT_EQ(descriptor.getDefaultOption(), "New dummy option");

  ASSERT_THROW(descriptor.setDefaultOption("Nonexisting option"), UniversalSettings::OptionDoesNotExistException);

  descriptor.setDefaultOption("Second new option");
  ASSERT_EQ(descriptor.getDefaultIndex(), 1);
  ASSERT_EQ(descriptor.getDefaultOption(), "Second new option");
}
} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
