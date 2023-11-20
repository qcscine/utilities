/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Typenames.h"
#include "gmock/gmock.h"

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {
/**
 * @class PositionCollectionTest PositionCollectionTest.cpp
 * @brief Comprises tests for the class Scine::Utils::PositionCollection.
 * @test
 */
class PositionCollectionTest : public Test {
 public:
  Position p1, p2, p3, p4;
  Displacement d1, d2, d3, d4;
  PositionCollection pc;
  DisplacementCollection dc;

  static void assertPositionCollectionsEqual(const PositionCollection& p1, const PositionCollection& p2);

 protected:
  void SetUp() override {
    d1 = Displacement(0.1, .2, -3);
    d2 = Displacement(0.3, 0.3, 0.1);
    d3 = Displacement(-0.3, 1.2, 0.1);
    d4 = Displacement(3.2, -2.2, 1);
    p1 = Position(1, 2, 3);
    p2 = Position(0.1, 0.2, 0.3);
    p3 = Position(-0.5, 0.2, 0.3);
    p4 = Position(9.2, -2.2, 2);
    pc.resize(4, 3);
    dc.resize(4, 3);
    pc.row(0) = p1;
    pc.row(1) = p2;
    pc.row(2) = p3;
    pc.row(3) = p4;
    dc.row(0) = d1;
    dc.row(1) = d2;
    dc.row(2) = d3;
    dc.row(3) = d4;
  }
};

TEST_F(PositionCollectionTest, CanBeAddedWithADisplacementCollection) {
  pc += dc;
  PositionCollection pc2(4, 3);
  pc2.row(0) = p1 + d1;
  pc2.row(1) = p2 + d2;
  pc2.row(2) = p3 + d3;
  pc2.row(3) = p4 + d4;
  assertPositionCollectionsEqual(pc, pc2);
}

TEST_F(PositionCollectionTest, CanBeMultipliedWithAFactor) {
  double f = 2.3;
  pc *= f;

  PositionCollection pc2(4, 3);
  pc2.row(0) = p1 * f;
  pc2.row(1) = p2 * f;
  pc2.row(2) = p3 * f;
  pc2.row(3) = p4 * f;

  assertPositionCollectionsEqual(pc, pc2);
}

void PositionCollectionTest::assertPositionCollectionsEqual(const PositionCollection& p1, const PositionCollection& p2) {
  ASSERT_THAT(p1.size(), Eq(p2.size()));
  for (int i = 0; i < p1.rows(); ++i) {
    ASSERT_THAT(p1.row(i).x(), DoubleEq(p2.row(i).x()));
    ASSERT_THAT(p1.row(i).y(), DoubleEq(p2.row(i).y()));
    ASSERT_THAT(p1.row(i).z(), DoubleEq(p2.row(i).z()));
  }
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
