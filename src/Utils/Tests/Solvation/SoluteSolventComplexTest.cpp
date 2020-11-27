/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h"
#include "Utils/IO/ChemicalFileFormats/XyzStreamHandler.h"
#include "Utils/Solvation/RandomIndexGenerator.h"
#include "Utils/Solvation/SoluteSolventComplex.h"
#include <gmock/gmock.h>
#include <fstream>
#include <iostream>
#include <random>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

// helper function to assert that Position Collections are equal.
static void assertPositionCollectionsEqual(const PositionCollection& p1, const PositionCollection& p2);

// Copied function from position collection test
static void assertPositionCollectionsEqual(const PositionCollection& p1, const PositionCollection& p2) {
  ASSERT_THAT(p1.size(), Eq(p2.size()));
  for (int i = 0; i < p1.rows(); ++i) {
    ASSERT_THAT(p1.row(i).x(), DoubleNear(p2.row(i).x(), 1e-15));
    ASSERT_THAT(p1.row(i).y(), DoubleNear(p2.row(i).y(), 1e-15));
    ASSERT_THAT(p1.row(i).z(), DoubleNear(p2.row(i).z(), 1e-15));
  }
}

TEST(SoluteSolventComplexTest, DoesCheckDistanceWork) {
  Atom testAtom(ElementType::Xe);

  AtomCollection atoms1;
  atoms1.push_back(testAtom);

  testAtom.setPosition(Position(ElementInfo::vdwRadius(testAtom.getElementType()) * 3, 0, 0));
  AtomCollection atomsNotOverlapping;
  atomsNotOverlapping.push_back(testAtom);
  ASSERT_THAT(SoluteSolventComplex::checkDistances(atoms1, atomsNotOverlapping), true);

  testAtom.setPosition(Position(ElementInfo::vdwRadius(testAtom.getElementType()) * 1.9, 0, 0));
  AtomCollection atomsOverlapping;
  atomsOverlapping.push_back(testAtom);
  ASSERT_THAT(SoluteSolventComplex::checkDistances(atoms1, atomsOverlapping), false);

  testAtom.setPosition(Position(ElementInfo::vdwRadius(testAtom.getElementType()) * 2.0, 0, 0));
  AtomCollection atomsTouching;
  atomsTouching.push_back(testAtom);
  ASSERT_THAT(SoluteSolventComplex::checkDistances(atoms1, atomsTouching), true);
}

// stack two tetraeder on each other, chosen vertice facing each other
TEST(SoluteSolventComplexTest, DoesArrangeAndRotateAroundAxisWorkProperly) {
  AtomCollection Td(4);
  // edge length of tetraeder is 2 * sqrt(2)
  ElementTypeCollection TdEC(4);
  TdEC.at(0) = ElementType::H;
  TdEC.at(1) = ElementType::He;
  TdEC.at(2) = ElementType::Li;
  TdEC.at(3) = ElementType::Be;

  // build Td around origin
  PositionCollection TdPC0(4, 3);
  TdPC0.row(0) = Position(sqrt(2), 0, -1);
  TdPC0.row(1) = Position(-sqrt(2), 0, -1);
  TdPC0.row(2) = Position(0, sqrt(2), 1);
  TdPC0.row(3) = Position(0, -sqrt(2), 1);

  Td.setElements(TdEC);
  // shift position to see if it works outside of the origin
  Position shiftVector = Position(-66, 42, 24);
  auto TdPC0_shifted = Geometry::translatePositions(TdPC0, shiftVector);
  Td.setPositions(TdPC0_shifted);

  int vertice = 0;

  auto TdPC1 =
      SoluteSolventComplex::arrange(TdPC0_shifted.row(vertice), (TdPC0_shifted.row(vertice) - shiftVector).normalized(),
                                    TdPC0.row(vertice + 1), TdPC0.row(vertice + 1).normalized(), TdPC0, 0);

  Td.setPositions(TdPC1);

  double distanceTh = 4 * sqrt(2) * cos(Constants::pi / 2 - 109.4712206 / 2 * Constants::pi / 180);

  // check that the distances are correct after arrangement
  ASSERT_THAT((TdPC0_shifted.row(vertice) - TdPC1.row(vertice + 1)).norm(), 0.0);
  ASSERT_THAT((TdPC0_shifted.row(vertice + 1) - TdPC1.row(vertice)).norm(), DoubleNear(distanceTh, 1e-9));
  ASSERT_THAT((TdPC0_shifted.row(2) - TdPC1.row(2)).norm(), DoubleNear(distanceTh, 1e-9));
  ASSERT_THAT((TdPC0_shifted.row(3) - TdPC1.row(3)).norm(), DoubleNear(distanceTh, 1e-9));

  // rotate by one 2pi * 1 / 3
  Position rotAxis = (TdPC0_shifted.row(vertice) - shiftVector).normalized();
  auto TdPC1_Rot1 = Geometry::rotatePositions(TdPC1, rotAxis, 2 * Constants::pi * 1 / 3, TdPC0_shifted.row(vertice));

  Td.setPositions(TdPC1_Rot1);

  // check that the distances are correct after rotation
  ASSERT_THAT((TdPC0_shifted.row(vertice) - TdPC1_Rot1.row(vertice + 1)).norm(), 0.0);
  ASSERT_THAT((TdPC0_shifted.row(1) - TdPC1_Rot1.row(3)).norm(), DoubleNear(distanceTh, 1e-9));
  ASSERT_THAT((TdPC0_shifted.row(2) - TdPC1_Rot1.row(0)).norm(), DoubleNear(distanceTh, 1e-9));
  ASSERT_THAT((TdPC0_shifted.row(3) - TdPC1_Rot1.row(2)).norm(), DoubleNear(distanceTh, 1e-9));
}

TEST(SoluteSolventComplexTest, DoesAddAlgorithmWorkCorrectly) {
  PositionCollection methanalPC(4, 3);
  methanalPC.row(0) = Position(2, 0, 0);
  methanalPC.row(1) = Position(0, 0, 0);
  methanalPC.row(2) = Position(-1.8, 1.5, 0);
  methanalPC.row(3) = Position(-1.8, -1.5, 0);

  ElementTypeCollection methanalEC(4);
  methanalEC.at(0) = ElementType::O;
  methanalEC.at(1) = ElementType::C;
  methanalEC.at(2) = ElementType::H;
  methanalEC.at(3) = ElementType::H;

  AtomCollection methanal(methanalEC, methanalPC);

  AtomCollection orgMethanal = methanal;

  ElementTypeCollection waterEC(3);
  waterEC.at(0) = ElementType::O;
  waterEC.at(1) = ElementType::H;
  waterEC.at(2) = ElementType::H;

  PositionCollection waterPC0(3, 3);
  waterPC0.row(0) = Position(0, 0, 0);
  waterPC0.row(1) = Position(-1.8, 1.25, 0);
  waterPC0.row(2) = Position(-1.8, -1.25, 0);

  AtomCollection water(waterEC, waterPC0);

  MolecularSurface::SurfaceSite attackSite(Position(0, 0, 0), Position(-1, 0, 0));

  double distance = 2;
  // minimal allowed distance
  double touchingDist = ElementInfo::vdwRadius(ElementType::O) * 2 + 2;
  // Test for no successful addition
  bool result1 = SoluteSolventComplex::add(methanal, water, attackSite, attackSite, distance, distance + 0.1, 0.25);
  methanal = orgMethanal;

  ASSERT_FALSE(result1);

  // Test for successful addition
  bool result2 = SoluteSolventComplex::add(methanal, water, attackSite, attackSite, touchingDist, touchingDist + 0.1, 0.25);
  // Distance between two oxygen atoms
  double controlDist2 = (methanal.at(orgMethanal.size()).getPosition() - methanal.at(1).getPosition()).norm();

  ASSERT_TRUE(result2);
  ASSERT_THAT(methanal.size(), orgMethanal.size() + water.size());
  ASSERT_THAT(controlDist2, touchingDist);
  methanal = orgMethanal;

  // Test for rotation of additive
  // add blocker
  methanal.push_back(Atom(ElementType::H, Position(12.1, -4.4, 0)));
  bool result3 = SoluteSolventComplex::add(methanal, water, attackSite, attackSite, touchingDist, touchingDist + 0.1, 0.25);
  double controlDist3 = (methanal.at(orgMethanal.size() + 1).getPosition() - methanal.at(1).getPosition()).norm();

  ASSERT_TRUE(result3);
  ASSERT_THAT(methanal.size(), orgMethanal.size() + 1 + water.size());
  ASSERT_THAT(controlDist3, touchingDist);
  ASSERT_THAT(methanal.at(6).getPosition().z(), DoubleNear(-1.0825, 1e-4));
}

TEST(SoluteSolventComplexTest, DoesStrategyDeciderWork) {
  ASSERT_THAT(SoluteSolventComplex::solvationStrategy(100000), 1000);
  ASSERT_THAT(SoluteSolventComplex::solvationStrategy(10000), 100);
  ASSERT_THAT(SoluteSolventComplex::solvationStrategy(501), 50);
  ASSERT_THAT(SoluteSolventComplex::solvationStrategy(499), 25);
  ASSERT_THAT(SoluteSolventComplex::solvationStrategy(201), 25);
  ASSERT_THAT(SoluteSolventComplex::solvationStrategy(200), 10);
  ASSERT_THAT(SoluteSolventComplex::solvationStrategy(101), 10);
  ASSERT_THAT(SoluteSolventComplex::solvationStrategy(100), 5);
  ASSERT_THAT(SoluteSolventComplex::solvationStrategy(21), 5);
  ASSERT_THAT(SoluteSolventComplex::solvationStrategy(20), 1);
  ASSERT_THAT(SoluteSolventComplex::solvationStrategy(0), 1);
}

// Test if Random Index Generator produces same number with same seed and if the number lies between 0 and size - 1
TEST(SoluteSolventComplexTest, DoesRandomIndexGeneratorWork) {
  int size = 100;
  int seed = 5;

  SoluteSolventComplex::RandomIndexGenerator generator1(size, seed);
  SoluteSolventComplex::RandomIndexGenerator generator2(size, seed);
  int i = 1;
  while (i < 1001) {
    int randInt1 = generator1.next();
    int randInt2 = generator2.next();

    ASSERT_THAT(randInt1, randInt2);
    ASSERT_THAT(true, randInt1 >= 0 && randInt1 < 100);
    i += 1;
  }
}

TEST(SoluteSolventComplexTest, DoesMergeAtomCollectionVectorWork) {
  AtomCollection test1;
  AtomCollection test2;

  test1.push_back(Atom(ElementType::He, Position(0, 0, 0)));
  test2.push_back(Atom(ElementType::Ar, Position(10, 10, 10)));

  int assembleSize = 10;
  std::vector<AtomCollection> assemble;

  for (int i = 0; i < assembleSize; i++) {
    if (i % 2 == 0) {
      assemble.push_back(test1);
    }
    else {
      assemble.push_back(test2);
    }
  }

  auto mergedAssemble = SoluteSolventComplex::mergeAtomCollectionVector(assemble);

  ASSERT_THAT(mergedAssemble.size(), assembleSize);
  ASSERT_TRUE(mergedAssemble.at(0).getElementType() == test1.at(0).getElementType());
  ASSERT_TRUE(mergedAssemble.at(mergedAssemble.size() - 1).getElementType() == test2.at(0).getElementType());
}

TEST(SoluteSolventComplexTest, DoesMergeSolventVectorWork) {
  AtomCollection test1;
  AtomCollection test2;

  test1.push_back(Atom(ElementType::He, Position(0, 0, 0)));
  test2.push_back(Atom(ElementType::Ar, Position(10, 10, 10)));

  int initAssembleSize = 1;
  int assembleCollSize = 5;

  std::vector<std::vector<AtomCollection>> assembleColl;

  for (int j = 0; j < assembleCollSize; j++) {
    std::vector<AtomCollection> assemble;
    for (int i = 0; i < initAssembleSize; i++) {
      if (i % 2 == 0) {
        assemble.push_back(test1);
      }
      else {
        assemble.push_back(test2);
      }
    }
    assembleColl.push_back(assemble);
    initAssembleSize *= 2;
  }

  auto mergedAssembleColl = SoluteSolventComplex::mergeSolventShellVector(assembleColl);

  ASSERT_THAT(mergedAssembleColl.size() + 1, pow(2, assembleCollSize));
  ASSERT_TRUE(mergedAssembleColl.at(0).getElementType() == test1.at(0).getElementType());
  // even indices should have test2, odd indices should have test1
  for (int i = 1; i < mergedAssembleColl.size() / 2; i++) {
    ASSERT_TRUE(mergedAssembleColl.at(2 * i - 1).getElementType() == test1.at(0).getElementType());
    ASSERT_TRUE(mergedAssembleColl.at(2 * i).getElementType() == test2.at(0).getElementType());
  }
}

TEST(SoluteSolventComplexTest, DoesSolvateAddsCorrectNumberOfSolvents) {
  ElementTypeCollection waterEC(3);
  waterEC.at(0) = ElementType::O;
  waterEC.at(1) = ElementType::H;
  waterEC.at(2) = ElementType::H;

  PositionCollection waterPC0(3, 3);
  waterPC0.row(0) = Position(0, 0, 0);
  waterPC0.row(1) = Position(-1.63, 1.25, 0);
  waterPC0.row(2) = Position(1.63, 1.25, 0);

  AtomCollection water(waterEC, waterPC0);

  AtomCollection solute;
  solute.push_back(Atom(ElementType::Na, Position(0, 0, 0)));

  int sampleRes = 18;
  int numSolvents = 24;

  auto solventComplex =
      SoluteSolventComplex::solvate(solute, solute.size(), water, numSolvents, 5, sampleRes, 0, 7, 1, 3, false);

  ASSERT_THAT(numSolvents, SoluteSolventComplex::mergeSolventShellVector(solventComplex).size() / water.size());

  // test for iteratively building up complex, should add 18 solvent molecules
  int i = 0;
  auto iterComplex = solute;

  while (i < 18) {
    auto newComplex = SoluteSolventComplex::solvate(iterComplex, solute.size(), water, 1, i + 5, sampleRes, 0, 7, 1, 3);
    iterComplex += SoluteSolventComplex::mergeSolventShellVector(newComplex);
    i++;
  }

  auto iterComplexSurface = MolecularSurface::getVisibleMolecularSurface(iterComplex, 0, solute.size(), sampleRes);

  ASSERT_THAT((iterComplex.size() - solute.size()) / water.size(), 18);
  ASSERT_THAT(iterComplexSurface.size(), 0);
}

TEST(SoluteSolventComplexTest, DoesSolvateFinishesSolventShell) {
  ElementTypeCollection waterEC(3);
  waterEC.at(0) = ElementType::O;
  waterEC.at(1) = ElementType::H;
  waterEC.at(2) = ElementType::H;

  PositionCollection waterPC0(3, 3);
  waterPC0.row(0) = Position(0, 0, 0);
  waterPC0.row(1) = Position(-1.63, 1.25, 0);
  waterPC0.row(2) = Position(1.63, 1.25, 0);

  AtomCollection water(waterEC, waterPC0);

  AtomCollection solute;
  solute.push_back(Atom(ElementType::Na, Position(0, 0, 0)));

  int sampleRes = 18;
  int numShells = 2;

  auto solventComplex =
      SoluteSolventComplex::solvateShells(solute, solute.size(), water, numShells, 5, sampleRes, 0, 7, 1, 3, false);
  auto soluteComplex = solute + SoluteSolventComplex::mergeSolventShellVector(solventComplex);

  auto soluteComplexSurface =
      MolecularSurface::getVisibleMolecularSurface(soluteComplex, solute.size(), solute.size() + 18 * 3, sampleRes);

  // first solvent shell is full after addition of 18 water molecules
  ASSERT_THAT(soluteComplexSurface.size(), 0);
  // first and second shell together contain approximatly 130 water molecules (Compiler dependent)
  // seed: 5, resolution: 18, solventOffset: 0, maxDistance: 7, stepSize: 1, numRotamers: 3, strategy: false; release mode
  ASSERT_THAT(true, (soluteComplex.size() - solute.size()) / water.size() > 2 * 18);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
