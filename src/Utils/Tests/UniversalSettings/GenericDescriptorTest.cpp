/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/UniversalSettings/BoolDescriptor.h>
#include <Utils/UniversalSettings/CollectionListDescriptor.h>
#include <Utils/UniversalSettings/DescriptorCollection.h>
#include <Utils/UniversalSettings/DirectoryDescriptor.h>
#include <Utils/UniversalSettings/DoubleDescriptor.h>
#include <Utils/UniversalSettings/DoubleListDescriptor.h>
#include <Utils/UniversalSettings/FileDescriptor.h>
#include <Utils/UniversalSettings/GenericDescriptor.h>
#include <Utils/UniversalSettings/IntDescriptor.h>
#include <Utils/UniversalSettings/IntListDescriptor.h>
#include <Utils/UniversalSettings/OptionListDescriptor.h>
#include <Utils/UniversalSettings/ParametrizedOptionListDescriptor.h>
#include <Utils/UniversalSettings/StringDescriptor.h>
#include <Utils/UniversalSettings/StringListDescriptor.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class GenericDescriptorTest GenericDescriptorTest.cpp
 * @brief Comprises tests for the class Scine::Utils::UniversalSettings::GenericDescriptor.
 * @test
 */
class GenericDescriptorTest : public Test {};

TEST_F(GenericDescriptorTest, ReturnsCorrectType) {
  auto boolDescriptor = UniversalSettings::BoolDescriptor("BoolDescriptor");
  auto descriptor = UniversalSettings::GenericDescriptor(boolDescriptor);
  ASSERT_TRUE(descriptor.relatesToBool());
  ASSERT_EQ(descriptor.getType(), UniversalSettings::GenericDescriptor::Type::Bool);

  auto intDescriptor = UniversalSettings::IntDescriptor("IntDescriptor");
  descriptor = UniversalSettings::GenericDescriptor(intDescriptor);
  ASSERT_TRUE(descriptor.relatesToInt());
  ASSERT_EQ(descriptor.getType(), UniversalSettings::GenericDescriptor::Type::Int);

  auto doubleDescriptor = UniversalSettings::DoubleDescriptor("DoubleDescriptor");
  descriptor = UniversalSettings::GenericDescriptor(doubleDescriptor);
  ASSERT_TRUE(descriptor.relatesToDouble());
  ASSERT_EQ(descriptor.getType(), UniversalSettings::GenericDescriptor::Type::Double);

  auto stringDescriptor = UniversalSettings::StringDescriptor("StringDescriptor");
  descriptor = UniversalSettings::GenericDescriptor(stringDescriptor);
  ASSERT_TRUE(descriptor.relatesToString());
  ASSERT_EQ(descriptor.getType(), UniversalSettings::GenericDescriptor::Type::String);

  auto fileDescriptor = UniversalSettings::FileDescriptor("FileDescriptor");
  descriptor = UniversalSettings::GenericDescriptor(fileDescriptor);
  ASSERT_TRUE(descriptor.relatesToFile());
  ASSERT_EQ(descriptor.getType(), UniversalSettings::GenericDescriptor::Type::File);

  auto directoryDescriptor = UniversalSettings::DirectoryDescriptor("DirectoryDescriptor");
  descriptor = UniversalSettings::GenericDescriptor(directoryDescriptor);
  ASSERT_TRUE(descriptor.relatesToDirectory());
  ASSERT_EQ(descriptor.getType(), UniversalSettings::GenericDescriptor::Type::Directory);

  auto optionListDescriptor = UniversalSettings::OptionListDescriptor("OptionListDescriptor");
  descriptor = UniversalSettings::GenericDescriptor(optionListDescriptor);
  ASSERT_TRUE(descriptor.relatesToOptionList());
  ASSERT_EQ(descriptor.getType(), UniversalSettings::GenericDescriptor::Type::OptionList);

  auto settingCollectionDescriptor = UniversalSettings::DescriptorCollection("SettingCollectionDescriptor");
  descriptor = UniversalSettings::GenericDescriptor(settingCollectionDescriptor);
  ASSERT_TRUE(descriptor.relatesToSettingCollection());
  ASSERT_EQ(descriptor.getType(), UniversalSettings::GenericDescriptor::Type::SettingCollection);

  auto parametrizedOptionListDescriptor =
      UniversalSettings::ParametrizedOptionListDescriptor("ParametrizedOptionListDescriptor");
  descriptor = UniversalSettings::GenericDescriptor(parametrizedOptionListDescriptor);
  ASSERT_TRUE(descriptor.relatesToParametrizedOptionList());
  ASSERT_EQ(descriptor.getType(), UniversalSettings::GenericDescriptor::Type::ParametrizedOptionList);

  auto intListDescriptor = UniversalSettings::IntListDescriptor("IntListDescriptor");
  descriptor = UniversalSettings::GenericDescriptor(intListDescriptor);
  ASSERT_TRUE(descriptor.relatesToIntList());
  ASSERT_EQ(descriptor.getType(), UniversalSettings::GenericDescriptor::Type::IntList);

  auto doubleListDescriptor = UniversalSettings::DoubleListDescriptor("DoubleListDescriptor");
  descriptor = UniversalSettings::GenericDescriptor(doubleListDescriptor);
  ASSERT_TRUE(descriptor.relatesToDoubleList());
  ASSERT_EQ(descriptor.getType(), UniversalSettings::GenericDescriptor::Type::DoubleList);

  auto stringListDescriptor = UniversalSettings::StringListDescriptor("StringListDescriptor");
  descriptor = UniversalSettings::GenericDescriptor(stringListDescriptor);
  ASSERT_TRUE(descriptor.relatesToStringList());
  ASSERT_EQ(descriptor.getType(), UniversalSettings::GenericDescriptor::Type::StringList);

  auto descriptorCollection = UniversalSettings::DescriptorCollection("DescriptorCollection");
  auto collectionListDescriptor =
      UniversalSettings::CollectionListDescriptor("CollectionListDescriptor", descriptorCollection);
  descriptor = UniversalSettings::GenericDescriptor(collectionListDescriptor);
  ASSERT_TRUE(descriptor.relatesToCollectionList());
  ASSERT_EQ(descriptor.getType(), UniversalSettings::GenericDescriptor::Type::CollectionList);
}
} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
