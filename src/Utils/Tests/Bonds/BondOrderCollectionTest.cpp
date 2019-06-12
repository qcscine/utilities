/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Bonds/BondOrderCollection.h"
#include <gmock/gmock.h>
#include <Eigen/SparseCore>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

TEST(BondOrderCollectionTest, Empty) {
  int size = 3;
  BondOrderCollection b;
  b.resize(3);
  ASSERT_EQ(b.getSystemSize(), size);
  ASSERT_EQ(b.getOrder(0, 0), 0);
  ASSERT_EQ(b.getOrder(0, 1), 0);
}

TEST(BondOrderCollectionTest, SetZero) {
  int size = 3;
  BondOrderCollection b;
  b.resize(3);
  b.setOrder(0, 1, 3.0);
  ASSERT_EQ(b.getOrder(0, 1), 3.0);
  ASSERT_EQ(b.getOrder(1, 0), 3.0);
  b.setZero();
  ASSERT_EQ(b.getOrder(0, 1), 0);
  ASSERT_EQ(b.getOrder(1, 0), 0);
}

TEST(BondOrderCollectionTest, SetMatrixAndCompare) {
  int size = 3;
  BondOrderCollection b1;
  BondOrderCollection b2;
  BondOrderCollection b3;
  Eigen::SparseMatrix<double> same(3, 3);
  Eigen::SparseMatrix<double> diff(3, 3);

  same.insert(0, 1) = 5.0;
  same.insert(1, 0) = 5.0;
  same.insert(2, 2) = 2.2;

  diff.insert(1, 1) = 4.0;
  diff.insert(2, 1) = 3.0;
  diff.insert(1, 2) = 3.0;

  b1.setMatrix(same);
  b2.setMatrix(same);
  b3.setMatrix(diff);

  ASSERT_TRUE(b1 == b2);
  ASSERT_TRUE(b1 != b3);
  ASSERT_FALSE(b1 == b3);
  ASSERT_FALSE(b1 != b2);
}

TEST(BondOrderCollectionTest, SizeOne) {
  BondOrderCollection b;
  b.resize(1);
  b.setOrder(0, 0, 1);
  ASSERT_EQ(b.getOrder(0, 0), 1);
}

TEST(BondOrderCollectionTest, SizeTwo) {
  BondOrderCollection b;
  b.resize(2);
  b.setOrder(0, 1, 1);
  ASSERT_EQ(b.getOrder(0, 1), 1);
  ASSERT_EQ(b.getOrder(1, 0), 1);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
