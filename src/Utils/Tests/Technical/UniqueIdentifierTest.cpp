/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/MSVCCompatibility.h"
#include <Utils/Technical/UniqueIdentifier.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

TEST(UniqueIdentifierTest, IsDifferentForEachInstance) {
  UniqueIdentifier id1;
  UniqueIdentifier id2;

  ASSERT_THAT(id1, Ne(id2));
}

TEST(UniqueIdentifierTest, CanBeCopied) {
  UniqueIdentifier id1;
  UniqueIdentifier id2 = id1;

  ASSERT_THAT(id1, Eq(id2));
}

TEST(UniqueIdentifierTest, CanBeMoved) {
  UniqueIdentifier original;
  UniqueIdentifier originalCopy = original;
  UniqueIdentifier moveCopy = std::move(original);

  ASSERT_THAT(moveCopy, Eq(originalCopy));
}

TEST(UniqueIdentifierTest, IsDifferentAfterHavingBeenMovedFrom) {
  UniqueIdentifier id1;
  UniqueIdentifier id2 = std::move(id1);

  ASSERT_THAT(id1, Ne(id2));
}

TEST(UniqueIdentifierTest, HasAStringRepresentation) {
  UniqueIdentifier id1;
  UniqueIdentifier id2;

  ASSERT_FALSE(id1.getStringRepresentation().empty());
  ASSERT_THAT(id1.getStringRepresentation(), Ne(id2.getStringRepresentation()));
}

TEST(UniqueIdentifierTest, AssignmentOperators) {
  UniqueIdentifier id1;
  UniqueIdentifier id2;

  id1 = id2;
  ASSERT_TRUE(id1 == id2);
  id1 = UniqueIdentifier();
  ASSERT_TRUE(id1 != id2);
}

TEST(UniqueIdentifierTest, ComparisonOperators) {
  UniqueIdentifier id1;
  UniqueIdentifier id2;

  ASSERT_FALSE(id1 == id2);
  ASSERT_TRUE(id1 != id2);
  ASSERT_TRUE((id1 < id2) xor (id2 < id1));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
