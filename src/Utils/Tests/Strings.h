/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Strings.h"
#include "gmock/gmock.h"

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class StringTest StringTest.cpp
 * @brief Comprises tests for the functions in Strings.h
 * @test
 */
class StringTest : public Test {
 public:
  std::string a, b;

 protected:
  void SetUp() override {
  }
};

TEST_F(StringTest, CaseInsensitiveWorks) {
  a = "test";
  b = "test";
  ASSERT_TRUE(caseInsensitiveEqual(a, b));
  b = "tEst";
  ASSERT_TRUE(caseInsensitiveEqual(a, b));
  b = "something";
  ASSERT_FALSE(caseInsensitiveEqual(a, b));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
