/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Typenames.h"
#include "gmock/gmock.h"

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

template<class VC1, class VC2>
bool identical(const VC1& vc1, const VC2& vc2) {
  if (vc1.size() != vc2.size())
    return false;
  for (int i = 0; i < vc1.rows(); ++i) {
    auto v1 = vc1.row(i);
    auto v2 = vc2.row(i);
    if (v1.x() != v2.x() || v1.y() != v2.y() || v1.z() != v2.z())
      return false;
  }
  return true;
}

/**
 * @class VectorCollectionTest VectorCollectionTest.cpp
 * @brief Comprises tests for the classes in Typenames.h .
 * @test
 */
class VectorCollectionTest : public Test {
 public:
  Position p1, p2, p3, p4;
  Gradient g1, g2, g3, g4;
  PositionCollection pc;
  GradientCollection gc;

 protected:
  void SetUp() override {
    p1 = Position(1, 2, 3);
    p2 = Position(0.1, 0.2, 0.3);
    p3 = Position(-0.5, 0.2, 0.3);
    p4 = Position(9.2, -2.2, 2);
    g1 = Gradient(0.1, .2, -3);
    g2 = Gradient(0.3, 0.3, 0.1);
    g3 = Gradient(-0.3, 1.2, 0.1);
    g4 = Gradient(3.2, -2.2, 1);
    pc.resize(4, 3);
    gc.resize(4, 3);
    pc.row(0) = p1;
    pc.row(1) = p2;
    pc.row(2) = p3;
    pc.row(3) = p4;
    gc.row(0) = g1;
    gc.row(1) = g2;
    gc.row(2) = g3;
    gc.row(3) = g4;
  }
};

TEST_F(VectorCollectionTest, CanBeAssignedFromAnotherOne) {
  GradientCollection assignedGc;
  assignedGc = pc;
  ASSERT_TRUE(identical(pc, assignedGc));
}

TEST_F(VectorCollectionTest, CanAddAnotherVectorCollection) {
  auto expected = gc;
  expected.row(0) += p1;
  expected.row(1) += p2;
  expected.row(2) += p3;
  expected.row(3) += p4;

  gc += pc;
  ASSERT_TRUE(identical(gc, expected));
}

TEST_F(VectorCollectionTest, CanSubtract) {
  auto expected = gc;
  expected.row(0) -= p1;
  expected.row(1) -= p2;
  expected.row(2) -= p3;
  expected.row(3) -= p4;

  gc -= pc;
  ASSERT_TRUE(identical(gc, expected));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
