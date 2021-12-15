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
