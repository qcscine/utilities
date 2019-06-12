/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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
  void TestBody() override {
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

  ASSERT_THAT(atoms.size(), Eq(randomElements.size()));
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

TEST_F(AtomCollectionTest, AllowsForRangeBasedLoops) {
  AtomCollection atoms(randomElements, randomPositions);

  int i = 0;
  for (const auto& a : atoms) {
    ASSERT_THAT(a.getElementType(), Eq(atoms[i].getElementType()));
    ASSERT_THAT(a.getPosition(), Eq(atoms[i].getPosition()));
    ++i;
  }
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
