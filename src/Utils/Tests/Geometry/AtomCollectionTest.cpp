/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry/AtomCollection.h"
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {
/**
 * @class AtomCollectionTest AtomCollectionTest.cpp
 * @brief Comprises tests for the class Scine::Utils::AtomCollection.
 * @test
 */
class AtomCollectionTest : public Test {
 public:
  ElementTypeCollection randomElements;
  PositionCollection randomPositions;

 private:
  void SetUp() final {
    randomElements = ElementTypeCollection{ElementType::Am, ElementType::Ce, ElementType::Ca};
    randomPositions = Eigen::MatrixX3d::Random(3, 3);
  }
};

TEST_F(AtomCollectionTest, HasSizeZeroByDefault) {
  AtomCollection atoms;
  ASSERT_THAT(atoms.size(), Eq(0));
}

TEST_F(AtomCollectionTest, HasSizeConstructor) {
  AtomCollection atoms(5);

  ASSERT_THAT(atoms.size(), Eq(5));
  for (auto&& atom : atoms) {
    ASSERT_THAT(atom.getElementType(), Eq(ElementType::none));
    ASSERT_TRUE(atom.getPosition().isZero());
  }
}

TEST_F(AtomCollectionTest, CanBeGeneratedFromElementTypeAndPositionCollections) {
  AtomCollection atoms(randomElements, randomPositions);

  ASSERT_TRUE(randomElements.size() == 3);
  ASSERT_TRUE(atoms.size() == randomElements.size());
  for (int i = 0; i < atoms.size(); ++i) {
    auto atom = atoms[i];
    ASSERT_THAT(atom.getElementType(), Eq(randomElements[i]));
    ASSERT_THAT(atom.getPosition(), Eq(randomPositions.row(i)));
  }
}

TEST_F(AtomCollectionTest, CanBeCleared) {
  AtomCollection atoms(randomElements, randomPositions);

  atoms.clear();
  ASSERT_THAT(atoms.size(), Eq(0));
}

TEST_F(AtomCollectionTest, CanBeResized) {
  AtomCollection atoms(randomElements, randomPositions);

  atoms.resize(2);
  ASSERT_THAT(atoms.size(), Eq(2));
  atoms.resize(10);
  ASSERT_THAT(atoms.size(), Eq(10));
}

TEST_F(AtomCollectionTest, AtomsCanBeAdded) {
  AtomCollection atoms(randomElements, randomPositions);
  ElementType randomElement = ElementType::Ds;
  Position randomPosition(-454.454, 0.3343, 2.222);

  Atom a(randomElement, randomPosition);

  atoms.push_back(a);
  ASSERT_THAT(atoms.size(), Eq(randomElements.size() + 1));
  ASSERT_THAT(atoms[atoms.size() - 1].getElementType(), Eq(randomElement));
  ASSERT_THAT(atoms[atoms.size() - 1].getPosition(), Eq(randomPosition));
}

TEST_F(AtomCollectionTest, AtomsCanBeSwapped) {
  AtomCollection atoms(randomElements, randomPositions);
  AtomCollection other(randomElements, randomPositions);
  ASSERT_TRUE(atoms.size() == 3);
  ASSERT_TRUE(other.size() == 3);
  ASSERT_TRUE(atoms == other);
  other.swapIndices(0, 1);
  ASSERT_FALSE(atoms == other);
  ASSERT_TRUE(atoms[0].getElementType() == other[1].getElementType());
  ASSERT_TRUE(atoms[1].getElementType() == other[0].getElementType());
  ASSERT_TRUE(atoms[0].getPosition().isApprox(other[1].getPosition()));
  ASSERT_TRUE(atoms[1].getPosition().isApprox(other[0].getPosition()));
}

TEST_F(AtomCollectionTest, AllowsForRangeBasedLoops) {
  AtomCollection atoms(randomElements, randomPositions);

  int i = 0;
  for (const auto& a : atoms) {
    ASSERT_THAT(a.getElementType(), Eq(atoms[i].getElementType()));
    ASSERT_THAT(a.getPosition(), Eq(atoms[i].getPosition()));
    ++i;
  }
  ASSERT_EQ(atoms.size(), i);
}

TEST_F(AtomCollectionTest, AtomCollectionsCanBeCombined) {
  AtomCollection first;
  first.push_back(Atom(ElementType::He, Position(0, 0, 0)));

  AtomCollection second;
  second.push_back(Atom(ElementType::Xe, Position(1, 1, 1)));

  AtomCollection merge = first + second;

  ASSERT_THAT(merge.size(), first.size() + second.size());
  ASSERT_THAT(merge.getElement(0), first.getElement(0));
  ASSERT_THAT(merge.getElement(merge.size() - 1), second.getElement(second.size() - 1));

  merge += first + second;

  ASSERT_THAT(merge.size(), 2 * first.size() + 2 * second.size());
  ASSERT_THAT(merge.getElement(0), first.getElement(0));
  ASSERT_THAT(merge.getElement(merge.size() - 1), second.getElement(second.size() - 1));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
