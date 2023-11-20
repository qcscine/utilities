/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/UniversalSettings/GenericValue.h"
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/OptionListDescriptor.h>
#include <Utils/UniversalSettings/ParametrizedOptionListDescriptor.h>
#include <Utils/UniversalSettings/ParametrizedOptionValue.h>
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Utils {
namespace Tests {

TEST(Settings, InvalidValueExceptions) {
  using namespace UniversalSettings;

  DescriptorCollection descriptors;
  IntDescriptor multiplicity("spin multiplicity");
  multiplicity.setMinimum(1);
  descriptors.push_back("multiplicity", multiplicity);

  ValueCollection values;
  values.addInt("multiplicity", 0);

  Settings settings{values, descriptors};

  ASSERT_FALSE(settings.valid());
  ASSERT_THROW(settings.throwIncorrectSettings(), InvalidSettingsException);
}

TEST(Settings, DivergingKeysAndEquality) {
  using namespace UniversalSettings;

  DescriptorCollection descriptors;
  IntDescriptor multiplicity("spin multiplicity");
  IntDescriptor charge("charge");
  StringDescriptor method("method");

  descriptors.push_back("multiplicity", multiplicity);
  descriptors.push_back("charge", charge);
  descriptors.push_back("method", method);

  ValueCollection values;
  values.addInt("multiplicity", 1);
  values.addInt("charge", 0);
  values.addString("method", "bla");

  ValueCollection otherValues;
  otherValues.addInt("multiplicity", 3);
  otherValues.addInt("charge", 0);
  otherValues.addString("method", "blubb");

  ValueCollection duplicateValues;
  duplicateValues.addInt("multiplicity", 1);
  duplicateValues.addInt("charge", 0);
  duplicateValues.addString("method", "bla");

  Settings settings{values, descriptors};
  Settings otherSettings{otherValues, descriptors};
  Settings duplicateSettings{duplicateValues, descriptors};

  ASSERT_TRUE(settings.valid());
  ASSERT_TRUE(otherSettings.valid());

  ASSERT_FALSE(values == otherValues);
  ASSERT_FALSE(duplicateValues == otherValues);
  ASSERT_TRUE(values == duplicateValues);

  auto otherDiverging = values.getDivergingKeys(otherValues);
  ASSERT_TRUE(otherDiverging.size() == 2);
  ASSERT_TRUE(values.getDivergingKeys(otherValues, true).size() == 1);
  ASSERT_TRUE(std::find(otherDiverging.begin(), otherDiverging.end(), "multiplicity") != otherDiverging.end());
  ASSERT_TRUE(std::find(otherDiverging.begin(), otherDiverging.end(), "method") != otherDiverging.end());
  ASSERT_TRUE(std::find(otherDiverging.begin(), otherDiverging.end(), "charge") == otherDiverging.end());

  ASSERT_TRUE(values.getDivergingKeys(duplicateValues).empty());
}

TEST(Settings, ItemAccess) {
  using namespace UniversalSettings;

  DescriptorCollection descriptors;
  const IntDescriptor multiplicity("spin multiplicity");
  const IntDescriptor charge("charge");
  const StringDescriptor method("method");

  descriptors.push_back("multiplicity", multiplicity);
  descriptors.push_back("charge", charge);
  descriptors.push_back("method", method);

  ValueCollection values;
  values.addInt("multiplicity", 1);
  values.addInt("charge", 0);
  values.addString("method", "bla");

  ValueCollection otherValues;
  otherValues.addInt("multiplicity", 3);
  otherValues.addInt("charge", 0);
  otherValues.addString("method", "blubb");

  values.addCollection("other", otherValues);

  ASSERT_TRUE(values.getCollection("other") == otherValues);

  const auto keys = values.getKeys();
  const auto items = values.items();
  std::vector<std::string> itemKeys;
  std::vector<GenericValue> itemValues;
  for (const auto& item : items) {
    itemKeys.push_back(item.first);
    itemValues.push_back(item.second);
  }
  ASSERT_TRUE(std::any_of(itemKeys.begin(), itemKeys.end(),
                          [&](const auto& key) { return std::find(keys.begin(), keys.end(), key) != keys.end(); }));
  ASSERT_TRUE(std::any_of(itemValues.begin(), itemValues.end(), [&](const auto& value) { return value == otherValues; }));

  const Settings otherSettings{otherValues, descriptors};
  descriptors.push_back("other", otherSettings.getDescriptorCollection());
  const Settings settings{values, descriptors};
  ASSERT_TRUE(settings.valid());
}

TEST(Settings, Merging) {
  using namespace UniversalSettings;

  DescriptorCollection descriptors;
  IntDescriptor multiplicity("spin multiplicity");
  IntDescriptor charge("charge");
  StringDescriptor method("method");

  descriptors.push_back("multiplicity", multiplicity);
  descriptors.push_back("charge", charge);
  descriptors.push_back("method", method);

  ValueCollection values;
  values.addInt("multiplicity", 1);
  values.addInt("charge", 0);
  values.addString("method", "bla");

  ValueCollection otherValues;
  otherValues.addInt("multiplicity", 3);
  otherValues.addInt("charge", 0);
  otherValues.addString("method", "blubb");

  Settings settings{values, descriptors};
  Settings otherSettings{otherValues, descriptors};

  ASSERT_TRUE(settings.valid());
  ASSERT_TRUE(otherSettings.valid());

  settings.merge(otherSettings);
  ASSERT_TRUE(settings.valueExists("multiplicity"));
  ASSERT_TRUE(settings.getInt("multiplicity") == otherSettings.getInt("multiplicity"));

  otherValues.addString("additional_setting", "something");
  Settings moreSettings{otherValues, descriptors};

  ASSERT_THROW(settings.merge(moreSettings), std::logic_error);
  settings.merge(moreSettings, true);
  ASSERT_TRUE(settings.getInt("multiplicity") == otherSettings.getInt("multiplicity"));
  ASSERT_FALSE(settings.valueExists("additional_setting"));

  settings.mergeAll(moreSettings);
  ASSERT_TRUE(settings.valueExists("additional_setting"));
  ASSERT_TRUE(settings.getString("additional_setting") == moreSettings.getString("additional_setting"));
}

TEST(Settings, CoerceStringCases) {
  using namespace UniversalSettings;

  DescriptorCollection descriptors;

  const std::string fruitKey = "fruit";
  OptionListDescriptor fruitDescriptor = [&]() {
    OptionListDescriptor descriptor;
    descriptor.addOption("apple");
    descriptor.addOption("banana");
    descriptor.addOption("cherry");
    return descriptor;
  }();
  descriptors.push_back(fruitKey, fruitDescriptor);

  const std::string shapeKey = "shape";
  ParametrizedOptionListDescriptor shapeDescriptor = [&]() {
    ParametrizedOptionListDescriptor descriptor;

    DescriptorCollection circleDescriptor;
    circleDescriptor.push_back("radius", DoubleDescriptor("radius"));
    descriptor.addOption("circle", circleDescriptor);

    DescriptorCollection rectangleDescriptor;
    rectangleDescriptor.push_back("x", DoubleDescriptor("x"));
    rectangleDescriptor.push_back("y", DoubleDescriptor("y"));
    descriptor.addOption("rectangle", rectangleDescriptor);

    return descriptor;
  }();
  descriptors.push_back(shapeKey, shapeDescriptor);

  ValueCollection values;
  values.addString(fruitKey, "Banana");

  ValueCollection circleValues;
  circleValues.addDouble("radius", 1.0);
  values.addOptionWithSettings(shapeKey, ParametrizedOptionValue{"CIRCLE", circleValues});

  Settings settings{values, descriptors};
  ASSERT_FALSE(settings.valid());
  settings.normalizeStringCases();
  ASSERT_EQ(settings.getString(fruitKey), "banana");
  ASSERT_EQ(settings.getOptionWithSettings(shapeKey).selectedOption, "circle");
  ASSERT_TRUE(settings.valid());
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
