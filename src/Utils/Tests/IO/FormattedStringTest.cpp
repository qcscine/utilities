/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/IO/FormattedString.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Utils {

using namespace testing;

TEST(FormattedStringTest, CorrectString) {
  auto string = format("%05d -- %s -- %+14.9f", 10, "foo", 19.123401);
  ASSERT_EQ(string, "00010 -- foo --  +19.123401000");
}

} // namespace Utils
} // namespace Scine
