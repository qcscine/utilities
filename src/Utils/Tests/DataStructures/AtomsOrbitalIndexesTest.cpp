/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Utils {

using namespace testing;

class AnAtomOrbitalIndexesTest : public Test {};

TEST_F(AnAtomOrbitalIndexesTest, CanCorrectlyBeConstructed) {
  auto indices0 = AtomsOrbitalsIndexes();
  auto indices10 = AtomsOrbitalsIndexes(10);
  ASSERT_EQ(indices0.getNAtoms(), 0);
  ASSERT_EQ(indices10.getNAtoms(), 10);
}

TEST_F(AnAtomOrbitalIndexesTest, ZeroSizeCanBeExtended) {
  auto indices = AtomsOrbitalsIndexes();
  ASSERT_EQ(indices.getNAtoms(), 0);
  ASSERT_EQ(indices.getNAtomicOrbitals(), 0);
  indices.addAtom(3);
  indices.addAtom(9);
  ASSERT_EQ(indices.getNAtoms(), 2);
  ASSERT_EQ(indices.getNAtomicOrbitals(), 12);
  ASSERT_EQ(indices.getNOrbitals(1), 9);
  ASSERT_EQ(indices.getFirstOrbitalIndex(1), 3);
}

TEST_F(AnAtomOrbitalIndexesTest, CanBeClearedAndResized) {
  auto indices = AtomsOrbitalsIndexes(4);
  ASSERT_EQ(indices.getNAtoms(), 4);
  indices.clear();
  ASSERT_EQ(indices.getNAtoms(), 0);
  indices.setSize(6);
  ASSERT_EQ(indices.getNAtoms(), 6);
  indices.addAtom(2);
  ASSERT_EQ(indices.getNAtoms(), 6);
  ASSERT_EQ(indices.getNAtomicOrbitals(), 2);
  ASSERT_EQ(indices.getNOrbitals(0), 2);
}

} // namespace Utils
} // namespace Scine