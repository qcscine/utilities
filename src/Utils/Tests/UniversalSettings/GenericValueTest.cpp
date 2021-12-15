/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
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
 * @brief Comprises tests for the class Scine::Utils::GenericValue.
 * @test
 */
class GenericValueTest : public Test {};

using namespace UniversalSettings;

TEST_F(GenericValueTest, ReturnsCorrectType) {
  bool boolValue = false;
  GenericValue gv = GenericValue::fromBool(boolValue);
  ASSERT_TRUE(gv.isBool());
  ASSERT_FALSE(gv.isInt());
  ASSERT_FALSE(gv.isDouble());
  ASSERT_FALSE(gv.isString());
  ASSERT_FALSE(gv.isCollection());
  ASSERT_FALSE(gv.isOptionWithSettings());

  int intValue = 0;
  gv = GenericValue::fromInt(intValue);
  ASSERT_FALSE(gv.isBool());
  ASSERT_TRUE(gv.isInt());
  ASSERT_FALSE(gv.isDouble());
  ASSERT_FALSE(gv.isString());
  ASSERT_FALSE(gv.isCollection());
  ASSERT_FALSE(gv.isOptionWithSettings());

  double doubleValue = 0.1;
  gv = GenericValue::fromDouble(doubleValue);
  ASSERT_FALSE(gv.isBool());
  ASSERT_FALSE(gv.isInt());
  ASSERT_TRUE(gv.isDouble());
  ASSERT_FALSE(gv.isString());
  ASSERT_FALSE(gv.isCollection());
  ASSERT_FALSE(gv.isOptionWithSettings());

  std::string stringValue = "test string";
  gv = GenericValue::fromString(stringValue);
  ASSERT_FALSE(gv.isBool());
  ASSERT_FALSE(gv.isInt());
  ASSERT_FALSE(gv.isDouble());
  ASSERT_TRUE(gv.isString());
  ASSERT_FALSE(gv.isCollection());
  ASSERT_FALSE(gv.isOptionWithSettings());

  ValueCollection valueCollection;
  gv = GenericValue::fromCollection(valueCollection);
  ASSERT_FALSE(gv.isBool());
  ASSERT_FALSE(gv.isInt());
  ASSERT_FALSE(gv.isDouble());
  ASSERT_FALSE(gv.isString());
  ASSERT_TRUE(gv.isCollection());
  ASSERT_FALSE(gv.isOptionWithSettings());

  ParametrizedOptionValue parametrizedOptionValue("dummy string", valueCollection);
  gv = GenericValue::fromOptionWithSettings(parametrizedOptionValue);
  ASSERT_FALSE(gv.isBool());
  ASSERT_FALSE(gv.isInt());
  ASSERT_FALSE(gv.isDouble());
  ASSERT_FALSE(gv.isString());
  ASSERT_FALSE(gv.isCollection());
  ASSERT_TRUE(gv.isOptionWithSettings());
}

TEST_F(GenericValueTest, IsComparable) {
  bool boolValue = false;
  GenericValue gv1 = GenericValue::fromBool(boolValue);
  GenericValue gv2 = GenericValue::fromBool(boolValue);
  boolValue = true;
  GenericValue gv3 = GenericValue::fromBool(boolValue);
  ASSERT_TRUE(gv1 == gv2);
  ASSERT_TRUE(gv1 != gv3);
  ASSERT_FALSE(gv2 == gv3);

  int intValue = 0;
  gv1 = GenericValue::fromInt(intValue);
  gv2 = GenericValue::fromInt(intValue);
  intValue = 10;
  gv3 = GenericValue::fromInt(intValue);
  ASSERT_TRUE(gv1 == gv2);
  ASSERT_TRUE(gv1 != gv3);
  ASSERT_FALSE(gv2 == gv3);

  double doubleValue = 0.1;
  gv1 = GenericValue::fromDouble(doubleValue);
  gv2 = GenericValue::fromDouble(doubleValue);
  doubleValue = 0.3;
  gv3 = GenericValue::fromDouble(doubleValue);
  ASSERT_TRUE(gv1 == gv2);
  ASSERT_TRUE(gv1 != gv3);
  ASSERT_FALSE(gv2 == gv3);

  std::string stringValue = "test string";
  gv1 = GenericValue::fromString(stringValue);
  gv2 = GenericValue::fromString(stringValue);
  stringValue = "other string";
  gv3 = GenericValue::fromString(stringValue);
  ASSERT_TRUE(gv1 == gv2);
  ASSERT_TRUE(gv1 != gv3);
  ASSERT_FALSE(gv2 == gv3);

  ValueCollection valueCollection;
  gv1 = GenericValue::fromCollection(valueCollection);
  gv2 = GenericValue::fromCollection(valueCollection);
  ASSERT_TRUE(gv1 == gv2);
}

TEST_F(GenericValueTest, MemberInterface) {
  GenericValue f{4};
  ASSERT_TRUE(f.isInt());
  ASSERT_TRUE(f.toInt() == 4);
  ASSERT_TRUE(f == 4);
  ASSERT_FALSE(f == 5);
  ASSERT_FALSE(f == 0.7);

  f = true;
  ASSERT_TRUE(f.isBool());
  ASSERT_TRUE(f == true);
  ASSERT_FALSE(f != true);
  ASSERT_FALSE(f == false);
  ASSERT_TRUE(f != false);
  ASSERT_FALSE(f == 4);

  f = std::vector<int>{{1, 2, 3}};
  ASSERT_TRUE(f.isIntList());
  ASSERT_TRUE(f.toIntList().at(1) == 2);

  // Implicit conversion and member comparison operators
  f = 4;
  const int x = f;
  ASSERT_EQ(f, x);

  // Expect throw when casting to non-selected member type
  double y = 0;
  ASSERT_THROW(y = f, std::runtime_error);
  ASSERT_EQ(y, 0);

  f = std::vector<double>{{1.0, 2.0, 3.0}};
  ASSERT_TRUE(f.isDoubleList());
  ASSERT_TRUE(f.toDoubleList().at(1) == 2.0);

  // const char helpers
  GenericValue g{"hello"};
  ASSERT_TRUE(g == "hello");
  g = "goodbye";
  ASSERT_TRUE(g == "goodbye");

  // Implicit constructor use for "variadic" function call
  auto fn = [](const GenericValue& /* v */) {};
  fn("Hi");
  fn(std::string{"Hi"});
  fn(4);
  fn(1.5);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
