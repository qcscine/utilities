/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/UniversalSettings/GenericValue.h>
#include <Utils/UniversalSettings/ParametrizedOptionValue.h>
#include <Utils/UniversalSettings/ValueCollection.h>
#include <gmock/gmock.h>
#include <string>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class GenericValueTest GenericValueTest.cpp
 * @brief Comprises tests for the class Scine::Utils::UniversalSettings::GenericValue.
 * @test
 */
class GenericValueTest : public Test {};

TEST_F(GenericValueTest, ReturnsCorrectType) {
  bool boolValue = false;
  UniversalSettings::GenericValue gv = UniversalSettings::GenericValue::fromBool(boolValue);
  ASSERT_TRUE(gv.isBool());
  ASSERT_FALSE(gv.isInt());
  ASSERT_FALSE(gv.isDouble());
  ASSERT_FALSE(gv.isString());
  ASSERT_FALSE(gv.isCollection());
  ASSERT_FALSE(gv.isOptionWithSettings());

  int intValue = 0;
  gv = UniversalSettings::GenericValue::fromInt(intValue);
  ASSERT_FALSE(gv.isBool());
  ASSERT_TRUE(gv.isInt());
  ASSERT_FALSE(gv.isDouble());
  ASSERT_FALSE(gv.isString());
  ASSERT_FALSE(gv.isCollection());
  ASSERT_FALSE(gv.isOptionWithSettings());

  double doubleValue = 0.1;
  gv = UniversalSettings::GenericValue::fromDouble(doubleValue);
  ASSERT_FALSE(gv.isBool());
  ASSERT_FALSE(gv.isInt());
  ASSERT_TRUE(gv.isDouble());
  ASSERT_FALSE(gv.isString());
  ASSERT_FALSE(gv.isCollection());
  ASSERT_FALSE(gv.isOptionWithSettings());

  std::string stringValue = "test string";
  gv = UniversalSettings::GenericValue::fromString(stringValue);
  ASSERT_FALSE(gv.isBool());
  ASSERT_FALSE(gv.isInt());
  ASSERT_FALSE(gv.isDouble());
  ASSERT_TRUE(gv.isString());
  ASSERT_FALSE(gv.isCollection());
  ASSERT_FALSE(gv.isOptionWithSettings());

  UniversalSettings::ValueCollection valueCollection;
  gv = UniversalSettings::GenericValue::fromCollection(valueCollection);
  ASSERT_FALSE(gv.isBool());
  ASSERT_FALSE(gv.isInt());
  ASSERT_FALSE(gv.isDouble());
  ASSERT_FALSE(gv.isString());
  ASSERT_TRUE(gv.isCollection());
  ASSERT_FALSE(gv.isOptionWithSettings());

  UniversalSettings::ParametrizedOptionValue parametrizedOptionValue("dummy string", valueCollection);
  gv = UniversalSettings::GenericValue::fromOptionWithSettings(parametrizedOptionValue);
  ASSERT_FALSE(gv.isBool());
  ASSERT_FALSE(gv.isInt());
  ASSERT_FALSE(gv.isDouble());
  ASSERT_FALSE(gv.isString());
  ASSERT_FALSE(gv.isCollection());
  ASSERT_TRUE(gv.isOptionWithSettings());
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
