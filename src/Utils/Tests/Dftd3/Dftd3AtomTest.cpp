/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Dftd3/Dftd3Atom.h"
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

class ADftd3AtomTest : public Test {
 public:
  Dftd3::Dftd3Atom atom;

  void SetUp() override {
  }
};

// Test that a coord. number can be correctly set for one atom.
TEST_F(ADftd3AtomTest, CanSetCoordinationNumber) {
  atom.setCoordinationNumber(1.25);
  EXPECT_THAT(atom.getCoordinationNumber(), DoubleNear(1.25, 1e-10));
}

// Test that an index can be correctly set for one atom.
TEST_F(ADftd3AtomTest, CanSetIndex) {
  atom.setIndex(1);
  ASSERT_THAT(atom.getIndex(), Eq(1));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
