/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
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
