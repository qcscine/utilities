/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry.h"
#include "Utils/Math/QuaternionFit.h"
#include "Utils/MolecularTrajectory.h"
#include "Utils/Typenames.h"
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <gmock/gmock.h>
#include <Eigen/Geometry>

using namespace testing;
namespace Scine {
namespace Utils {
using namespace Geometry;
namespace Tests {

/**
 * @class Scine::Utils::Tests::GeometryTest GeometryTest.cpp
 * @brief Comprises tests for the class Scine::Utils::Geometry.
 * @test
 */
TEST(GeometryTest, DoesShiftPositionCollectionWork) {
  PositionCollection testPC(2, 3);
  testPC.row(0) = Position(5, 7, 6);
  testPC.row(1) = Position(-1, -1, -1);

  Position shiftVector = Position(1.0, 1.0, 1.0);

  PositionCollection translatedTestPC;
  translatedTestPC = Manipulations::translatePositions(testPC, shiftVector);

  testPC.row(0) = Position(6, 8, 7);
  testPC.row(1) = Position(0, 0, 0);

  for (int i = 0; i < translatedTestPC.rows(); ++i) {
    ASSERT_THAT(translatedTestPC.row(i).x(), DoubleNear(testPC.row(i).x(), 1e-12));
    ASSERT_THAT(translatedTestPC.row(i).y(), DoubleNear(testPC.row(i).y(), 1e-12));
    ASSERT_THAT(translatedTestPC.row(i).z(), DoubleNear(testPC.row(i).z(), 1e-12));
  }
}

TEST(GeometryTest, CanTransformAPositionVectorToAMatrix) {
  int numberParticles = 11;
  int numberDimensions = 3 * numberParticles;

  Eigen::VectorXd positionVector = Eigen::VectorXd::Random(numberDimensions);

  auto matrix = Transformations::positionVectorToMatrix(positionVector);

  for (int i = 0; i < numberParticles; ++i) {
    for (int j = 0; j < 3; ++j) {
      ASSERT_THAT(matrix(i, j), DoubleEq(positionVector(3 * i + j)));
    }
  }
}

TEST(GeometryTest, CanTransformAPositionMatrixToAVector) {
  int numberParticles = 11;

  Eigen::MatrixXd positionMatrix = Eigen::MatrixXd::Random(numberParticles, 3);

  auto vector = Transformations::positionMatrixToVector(positionMatrix);

  for (int i = 0; i < numberParticles; ++i) {
    for (int j = 0; j < 3; ++j) {
      ASSERT_THAT(vector(3 * i + j), DoubleEq(positionMatrix(i, j)));
    }
  }
}

TEST(GeometryTest, DoesRotatePositionCollectionWork) {
  PositionCollection testPC(3, 3);
  testPC.row(0) = Position(2, 2, 1);
  testPC.row(1) = Position(2.0, -2.0, 1.0);
  testPC.row(2) = Position(2.0, 0, -1.0);

  Position startOrientation = Position(-1.0, 0.0, 0.0);
  Position endOrientation = Position(0.0, -1.0, 0.0);
  Position rotOrigin = Position(0.0, 0.0, 0.0);

  PositionCollection rotatedPC = Manipulations::rotatePositions(testPC, startOrientation, endOrientation, rotOrigin);
  // manually rotated position collection
  testPC.row(0) = Position(-2, 2, 1);
  testPC.row(1) = Position(2.0, 2.0, 1.0);
  testPC.row(2) = Position(0, 2, -1.0);

  for (int i = 0; i < rotatedPC.rows(); ++i) {
    ASSERT_THAT(rotatedPC.row(i).x(), DoubleNear(testPC.row(i).x(), 1e-12));
    ASSERT_THAT(rotatedPC.row(i).y(), DoubleNear(testPC.row(i).y(), 1e-12));
    ASSERT_THAT(rotatedPC.row(i).z(), DoubleNear(testPC.row(i).z(), 1e-12));
  }
}

TEST(GeometryTest, DoesRotateAroundAxisRotateRightAngle) {
  PositionCollection waterPC0(3, 3);
  waterPC0.row(0) = Position(0, 0, 0);
  waterPC0.row(1) = Position(-1.8, 1.25, 0);
  waterPC0.row(2) = Position(-1.8, -1.25, 0);

  // shift waterPC0 somewhere
  Position axisOrg1(42, 24, -2.4);
  Manipulations::translatePositionsInPlace(waterPC0, axisOrg1);

  // Rotation around z-axis of 120 degrees
  double rotAngle1 = 120 * Constants::pi / 180;
  Displacement axis1(0, 0, 1);

  PositionCollection waterPC1 = Manipulations::rotatePositions(waterPC0, axis1, rotAngle1, axisOrg1);
  Position v1 = waterPC0.row(1) - waterPC0.row(0);
  Position v2 = waterPC1.row(1) - waterPC1.row(0);
  double dotPr1 = v1.dot(v2);
  double measuredAngle1 = acos(dotPr1 / (v1.norm() * v2.norm())) * 180 / Constants::pi;

  // Rotation around x-axis of 90 degrees
  double rotAngle2 = 90 * Constants::pi / 180;
  Displacement axis2(1, 0, 0);

  PositionCollection waterPC2 = Manipulations::rotatePositions(waterPC0, axis2, rotAngle2, axisOrg1);
  Position v3 = waterPC2.row(1) - waterPC0.row(1);
  Position v4 = waterPC2.row(2) - waterPC0.row(1);
  double dotPr2 = v3.dot(v4);
  double measuredAngle2 = acos(dotPr2 / (v3.norm() * v4.norm())) * 180 / Constants::pi;

  ASSERT_THAT(measuredAngle1, DoubleNear(120.0, 1e-12));
  ASSERT_THAT(measuredAngle2, DoubleNear(90.0, 1e-12));
}

TEST(GeometryTest, CanRandomlyDisplaceInPlace) {
  PositionCollection pc1(2, 3);
  PositionCollection pc2(2, 3);
  pc1 << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
  pc2 << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
  double maxDisplacement = 0.1;
  Manipulations::randomDisplacementInPlace(pc1, maxDisplacement);
  for (int i = 0; i < pc1.rows(); ++i) {
    for (int j = 0; j < pc1.cols(); ++j) {
      ASSERT_LE(std::fabs(pc1(i, j) - pc2(i, j)), maxDisplacement);
    }
  }
}

TEST(GeometryTest, CanGenerateRandomTrajectoryFromAtomCollectionAndBeReproduced) {
  PositionCollection pc(2, 3);
  pc << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
  auto ec = ElementTypeCollection{ElementType::H, ElementType::H};
  AtomCollection atoms(ec, pc);
  int numFrames = 10;
  double maxDisplacement = 0.1;
  MolecularTrajectory trajectory1 = Manipulations::randomDisplacementTrajectory(atoms, numFrames, maxDisplacement, 42);
  EXPECT_THAT(trajectory1.molecularSize(), pc.rows());
  EXPECT_THAT(trajectory1.size(), numFrames);
  /* Check if maxDisplacement worked */
  for (int i = 0; i < trajectory1.size(); ++i) {
    for (int j = 0; j < pc.rows(); ++j) {
      for (int k = 0; k < pc.cols(); ++k) {
        EXPECT_LE(std::fabs(trajectory1[i](j, k) - pc(j, k)), maxDisplacement);
      }
    }
  }
  MolecularTrajectory trajectory2 = Manipulations::randomDisplacementTrajectory(atoms, numFrames, maxDisplacement, 42);
  EXPECT_THAT(trajectory2.molecularSize(), pc.rows());
  EXPECT_THAT(trajectory2.size(), numFrames);
  /* Check if trajectories are identical */
  for (int i = 0; i < trajectory1.size(); ++i) {
    for (int j = 0; j < pc.rows(); ++j) {
      for (int k = 0; k < pc.cols(); ++k) {
        EXPECT_LE(std::fabs(trajectory1[i](j, k) - trajectory2[i](j, k)), 1e-12);
      }
    }
  }
}

TEST(GeometryTest, CanAlignTwoPositionCollections) {
  int numberParticles = 11;
  PositionCollection referencePositions = Eigen::MatrixXd::Random(numberParticles, 3);
  PositionCollection positions = Eigen::MatrixXd::Random(numberParticles, 3);

  Manipulations::alignPositions(referencePositions, positions);

  QuaternionFit fit(referencePositions, positions);

  Eigen::Vector3d translationVector = fit.getTransVector();
  Eigen::Matrix3d rotation = fit.getRotationMatrix();

  Eigen::AngleAxisd angleRotation{rotation};

  ASSERT_THAT(angleRotation.angle(), DoubleNear(0.0, 1e-6));
  ASSERT_TRUE(translationVector.isZero(1e-6));
}

TEST(GeometryTest, GetsCorrectMassVector) {
  auto ec = ElementTypeCollection{ElementType::H, ElementType::C, ElementType::Fe};

  auto masses = Properties::getMasses(ec);
  std::vector<double> expectedMasses = {ElementInfo::mass(ElementType::H), ElementInfo::mass(ElementType::C),
                                        ElementInfo::mass(ElementType::Fe)};

  ASSERT_THAT(masses, Eq(expectedMasses));
}

TEST(GeometryTest, AtomIsCorrectlyFoundInStructure) {
  std::stringstream ss("5\n\n"
                       "C      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.6287000000    0.6287000000    0.6287000000\n"
                       "H     -0.6287000000   -0.6287000000    0.6287000000\n"
                       "H     -0.6287000000    0.6287000000   -0.6287000000\n"
                       "H      0.6287000000   -0.6287000000   -0.6287000000\n");

  auto structure = XyzStreamHandler::read(ss);

  Position p1(-0.6287, -0.6287, 0.6287);
  p1 *= Constants::bohr_per_angstrom;
  Atom atom1(ElementType::H, p1);

  auto index = Distances::getIndexOfAtomInStructure(structure, atom1);
  ASSERT_THAT(index, Eq(2));

  Position p2(-0.6287, -0.7287, 0.6287);
  p2 *= Constants::bohr_per_angstrom;
  Atom atom2(ElementType::H, p2);
  Atom atom3(ElementType::Cl, p1);

  EXPECT_THROW(Distances::getIndexOfAtomInStructure(structure, atom2), std::runtime_error);
  EXPECT_THROW(Distances::getIndexOfAtomInStructure(structure, atom3), std::runtime_error);

  Position p3(-0.6287, -0.63, 0.6287);
  p3 *= Constants::bohr_per_angstrom;
  Atom atom4(ElementType::H, p3);

  EXPECT_THROW(Distances::getIndexOfAtomInStructure(structure, atom3), std::runtime_error);
  ASSERT_THAT(Distances::getIndexOfAtomInStructure(structure, atom4, 0.1), Eq(2));
}

TEST(GeometryTest, ClosestAtomIsCorrectlyIdentified) {
  std::stringstream ss("5\n\n"
                       "C      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.6287000000    0.6287000000    0.6287000000\n"
                       "H     -0.6287000000   -0.6287000000    0.6287000000\n"
                       "H     -0.6287000000    0.6287000000   -0.6287000000\n"
                       "H      0.6287000000   -0.6287000000   -0.6287000000\n");

  auto structure = XyzStreamHandler::read(ss);
  const auto& positions = structure.getPositions();

  Position p1(-0.6287, -0.6287, 0.6287);
  p1 *= Constants::bohr_per_angstrom;

  ASSERT_THAT(Distances::getIndexOfClosestAtom(positions, p1), Eq(2));
  ASSERT_THAT(Distances::getIndexOfClosestAtom(positions, p1, -0.001), Eq(2));
  ASSERT_THAT(Distances::getIndexOfClosestAtom(positions, p1, 1e-4), Eq(0));

  Position p2(-1.0, 1.0, -1.0);
  p2 *= Constants::bohr_per_angstrom;

  ASSERT_THAT(Distances::getIndexOfClosestAtom(positions, p2), Eq(3));
  ASSERT_THAT(Distances::getIndexOfClosestAtom(positions, p2, -0.001), Eq(3));
  ASSERT_THAT(Distances::getIndexOfClosestAtom(positions, p2, 1e-4), Eq(3));
}

TEST(GeometryTest, CorrectlyReturnsIndicesCloseToAtom) {
  std::stringstream ss("8\n\n"
                       "H 0.0 0.0 0.0\n"
                       "H 1.0 0.0 0.0\n"
                       "H 2.0 0.0 0.0\n"
                       "H 3.0 0.0 0.0\n"
                       "H 4.0 0.0 0.0\n"
                       "H 5.0 0.0 0.0\n"
                       "H 6.0 0.0 0.0\n"
                       "H 7.0 0.0 0.0");
  auto structure = XyzStreamHandler::read(ss);
  const PositionCollection& positions = structure.getPositions();
  double distanceThreshold = 2.1 * Constants::bohr_per_angstrom;
  int targetAtom = 3;
  bool onlyUpperTriangle = true;
  bool includeCurrentAtom = false;
  auto closeIndices =
      Distances::getIndicesCloseToAtom(positions, targetAtom, distanceThreshold, includeCurrentAtom, onlyUpperTriangle);

  EXPECT_EQ(closeIndices.size(), 2);
  EXPECT_EQ(closeIndices[0], 4);
  EXPECT_EQ(closeIndices[1], 5);

  onlyUpperTriangle = true;
  includeCurrentAtom = true;
  closeIndices =
      Distances::getIndicesCloseToAtom(positions, targetAtom, distanceThreshold, includeCurrentAtom, onlyUpperTriangle);

  EXPECT_EQ(closeIndices.size(), 3);
  EXPECT_EQ(closeIndices[0], 3);
  EXPECT_EQ(closeIndices[1], 4);
  EXPECT_EQ(closeIndices[2], 5);

  onlyUpperTriangle = false;
  includeCurrentAtom = true;
  closeIndices =
      Distances::getIndicesCloseToAtom(positions, targetAtom, distanceThreshold, includeCurrentAtom, onlyUpperTriangle);

  EXPECT_EQ(closeIndices.size(), 5);
  EXPECT_EQ(closeIndices[0], 1);
  EXPECT_EQ(closeIndices[1], 2);
  EXPECT_EQ(closeIndices[2], 3);
  EXPECT_EQ(closeIndices[3], 4);
  EXPECT_EQ(closeIndices[4], 5);
}

TEST(GeometryTest, CorrectlyConstructAtomPairList) {
  std::stringstream ss("8\n\n"
                       "H 0.0 0.0 0.0\n"
                       "H 1.0 0.0 0.0\n"
                       "H 2.0 0.0 0.0\n"
                       "H 3.0 0.0 0.0\n"
                       "H 4.0 0.0 0.0\n"
                       "H 5.0 0.0 0.0\n"
                       "H 6.0 0.0 0.0\n"
                       "H 7.0 0.0 0.0");
  auto structure = XyzStreamHandler::read(ss);
  const PositionCollection& positions = structure.getPositions();
  double distanceThreshold = 2.1 * Constants::bohr_per_angstrom;
  bool onlyUpperTriangle = true;
  bool sameAtomIsIncluded = false;

  auto atomPairList = Distances::constructAtomPairList(positions, distanceThreshold, sameAtomIsIncluded, onlyUpperTriangle);

  EXPECT_EQ(atomPairList.size(), 8);
  for (const auto& atomList : atomPairList) {
    if (atomList.first == 6)
      EXPECT_EQ(atomList.second.size(), 1);
    else if (atomList.first == 7)
      EXPECT_EQ(atomList.second.size(), 0);
    else
      EXPECT_EQ(atomList.second.size(), 2);
  }
  EXPECT_EQ(atomPairList[0][0], 1);
  EXPECT_EQ(atomPairList[0][1], 2);
  EXPECT_EQ(atomPairList[1][0], 2);
  EXPECT_EQ(atomPairList[1][1], 3);
  EXPECT_EQ(atomPairList[2][0], 3);
  EXPECT_EQ(atomPairList[2][1], 4);
  EXPECT_EQ(atomPairList[3][0], 4);
  EXPECT_EQ(atomPairList[3][1], 5);
  EXPECT_EQ(atomPairList[4][0], 5);
  EXPECT_EQ(atomPairList[4][1], 6);
  EXPECT_EQ(atomPairList[5][0], 6);
  EXPECT_EQ(atomPairList[5][1], 7);
  EXPECT_EQ(atomPairList[6][0], 7);

  sameAtomIsIncluded = true;
  atomPairList = Distances::constructAtomPairList(positions, distanceThreshold, sameAtomIsIncluded, onlyUpperTriangle);
  EXPECT_EQ(atomPairList.size(), 8);
  for (const auto& atomList : atomPairList) {
    if (atomList.first == 6)
      EXPECT_EQ(atomList.second.size(), 2);
    else if (atomList.first == 7)
      EXPECT_EQ(atomList.second.size(), 1);
    else
      EXPECT_EQ(atomList.second.size(), 3);
  }
  EXPECT_EQ(atomPairList[0][0], 0);
  EXPECT_EQ(atomPairList[0][1], 1);
  EXPECT_EQ(atomPairList[0][2], 2);
  EXPECT_EQ(atomPairList[1][0], 1);
  EXPECT_EQ(atomPairList[1][1], 2);
  EXPECT_EQ(atomPairList[1][2], 3);
  EXPECT_EQ(atomPairList[2][0], 2);
  EXPECT_EQ(atomPairList[2][1], 3);
  EXPECT_EQ(atomPairList[2][2], 4);
  EXPECT_EQ(atomPairList[3][0], 3);
  EXPECT_EQ(atomPairList[3][1], 4);
  EXPECT_EQ(atomPairList[3][2], 5);
  EXPECT_EQ(atomPairList[4][0], 4);
  EXPECT_EQ(atomPairList[4][1], 5);
  EXPECT_EQ(atomPairList[4][2], 6);
  EXPECT_EQ(atomPairList[5][0], 5);
  EXPECT_EQ(atomPairList[5][1], 6);
  EXPECT_EQ(atomPairList[5][2], 7);
  EXPECT_EQ(atomPairList[6][0], 6);
  EXPECT_EQ(atomPairList[6][1], 7);
  EXPECT_EQ(atomPairList[7][0], 7);
}

TEST(GeometryTest, Centroid) {
  PositionCollection diametric(2, 3);
  diametric << 1.0, 2.0, 3.0, -1.0, 2.0, -3.0;
  ASSERT_TRUE(Properties::getAveragePosition(diametric).isApprox(Position(0.0, 2.0, 0.0), 1e-10));
}

TEST(GeometryTest, DistanceSquaredCorrectlyCalculated) {
  Position position1(0, 0, 0);
  Position position2(2, 1, 1);
  EXPECT_DOUBLE_EQ(Distances::distance(position1, position2), std::sqrt(6.0));
  EXPECT_DOUBLE_EQ(Distances::distanceSquared(position1, position2), 6.0);
}

TEST(GeometryTest, DistanceCorrectlyCalculated) {
  PositionCollection positions1(1, 3);
  PositionCollection positions2(1, 3);
  Position position1(0, 0, 0);
  Position position2(2, 1, 1);
  positions1.row(0) = position1;
  positions2.row(0) = position2;
  EXPECT_DOUBLE_EQ(Distances::distance(positions1, positions2), std::sqrt(6.0));
}

TEST(GeometryTest, CanRecoverListOfDivergingAtomPositions) {
  PositionCollection reference(5, 3);
  reference << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0, 4.0;
  PositionCollection testPositions = reference;
  EXPECT_EQ(Distances::getListOfDivergingAtoms(reference, testPositions, 0.11).size(), 0);

  DisplacementCollection deviation = DisplacementCollection::Zero(5, 3);
  deviation(3, 2) -= 0.5;
  deviation(4, 2) += 0.5;
  testPositions += deviation;

  auto indices = Distances::getListOfDivergingAtoms(reference, testPositions, 0.11);
  EXPECT_EQ(indices.size(), 2);
  EXPECT_EQ(indices[0], 3);
  EXPECT_EQ(indices[1], 4);

  indices = Distances::getListOfDivergingAtoms(
      reference, testPositions, 0.11, {ElementType::C, ElementType::C, ElementType::C, ElementType::C, ElementType::C});
  EXPECT_EQ(indices.size(), 2);
  EXPECT_EQ(indices[0], 3);
  EXPECT_EQ(indices[1], 4);
}

TEST(GeometryTest, CanRecoverListOfIterativeDivergingAtomPositions) {
  PositionCollection reference(5, 3);
  reference << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0, 4.0;
  PositionCollection testPositions = reference;
  EXPECT_EQ(Distances::getListOfDivergingAtomsRobust(reference, testPositions, 0.11, 1.0).size(), 0);

  DisplacementCollection deviation = DisplacementCollection::Zero(5, 3);
  deviation(2, 2) += 4.0;
  testPositions += deviation;

  auto indices = Distances::getListOfDivergingAtomsRobust(reference, testPositions, 0.05, 1e-2);
  EXPECT_EQ(indices.size(), 1);
  EXPECT_EQ(indices[0], 2);

  indices = Distances::getListOfDivergingAtomsRobust(
      reference, testPositions, 0.05, 1e-2, 20,
      {ElementType::C, ElementType::C, ElementType::C, ElementType::C, ElementType::H});
  EXPECT_EQ(indices.size(), 1);
  EXPECT_EQ(indices[0], 2);

  std::stringstream struct1("60\n\n"
                            "C     13.7720199616   -2.6605442099   -3.3601390949\n"
                            "C     13.4654423517   -1.4277613397   -4.2011238303\n"
                            "C     12.0760821738   -1.4844307593   -4.8468448268\n"
                            "C     11.7713880441   -0.2434583772   -5.6916659471\n"
                            "C     10.3791497395   -0.3022614777   -6.3370181015\n"
                            "C     10.0993997813    0.8943917103   -7.2070508979\n"
                            "C      8.9606917610    1.5989243626   -7.2664238042\n"
                            "C      7.7055283284    1.3495327804   -6.4764313115\n"
                            "C      7.0968673612    2.6447088122   -5.8879868658\n"
                            "H      7.8602302192    3.1476425416   -5.2474874513\n"
                            "O      5.9868788428    2.2567434799   -5.0683758300\n"
                            "C      6.6284447425    3.6150758523   -6.9533759995\n"
                            "C      7.2139005258    4.8043982031   -7.1987662781\n"
                            "C      6.7528019823    5.7242942271   -8.2312778176\n"
                            "C      7.2802369280    6.9311578748   -8.5224672451\n"
                            "C      8.4428428577    7.5981438897   -7.8355690078\n"
                            "C      9.4970081304    8.0584279448   -8.8149216231\n"
                            "C     10.1008704751    9.2541596859   -8.8416743615\n"
                            "C      9.8655419638   10.4126098563   -7.9088116983\n"
                            "C      9.7341105761   11.7458102252   -8.6592405887\n"
                            "C      9.5509675264   12.9355968554   -7.7015768954\n"
                            "C      9.3635185557   14.2358673430   -8.4411141363\n"
                            "O      8.2033481673   14.2662388777   -9.1652086274\n"
                            "C     10.1880719840   15.3076079531   -8.4482640052\n"
                            "C     11.4957357368   15.4283316784   -7.7173292588\n"
                            "H     14.7627326033   -2.5911972580   -2.9159393389\n"
                            "H     13.7377627474   -3.5646289275   -3.9648117118\n"
                            "H     13.0497880103   -2.7732275850   -2.5540333087\n"
                            "H     14.2256065376   -1.3199970132   -4.9822169317\n"
                            "H     13.5398272453   -0.5314248735   -3.5762054645\n"
                            "H     11.3152208903   -1.5921920279   -4.0662806017\n"
                            "H     12.0001546339   -2.3809514664   -5.4718626144\n"
                            "H     12.5314599717   -0.1352425156   -6.4734975748\n"
                            "H     11.8473392518    0.6523770580   -5.0664930250\n"
                            "H      9.6159226132   -0.4046097748   -5.5594330508\n"
                            "H     10.3123420094   -1.2135875568   -6.9486121951\n"
                            "H     10.9195939032    1.1961939168   -7.8597856707\n"
                            "H      8.9135733563    2.4353481456   -7.9650004249\n"
                            "H      7.8874162734    0.6503017927   -5.6560573924\n"
                            "H      6.9526160026    0.8778407145   -7.1220429910\n"
                            "H      5.5035424602    3.0477361581   -4.7861228378\n"
                            "H      5.7628013690    3.3133029789   -7.5439623565\n"
                            "H      8.0785889542    5.0985340917   -6.6027050126\n"
                            "H      5.8968325838    5.3862068045   -8.8159423811\n"
                            "H      6.8247205396    7.5045409118   -9.3296836194\n"
                            "H      8.8990845799    6.9128402534   -7.1098494171\n"
                            "H      8.0721142024    8.4538133602   -7.2553766683\n"
                            "H      9.7921982560    7.3199118138   -9.5600508884\n"
                            "H     10.8569030356    9.4249133801   -9.6092036500\n"
                            "H     10.7084416173   10.4837993774   -7.2062092737\n"
                            "H      8.9742020591   10.2455582405   -7.2958522199\n"
                            "H      8.8831749553   11.7027674393   -9.3444371241\n"
                            "H     10.6237623767   11.9070925258   -9.2771635774\n"
                            "H     10.4059448339   13.0000392473   -7.0235716132\n"
                            "H      8.6672922981   12.7557213907   -7.0773996733\n"
                            "H     11.7287562127   14.5345394671   -7.1448059125\n"
                            "H      8.1238612138   15.1140630870   -9.6272266876\n"
                            "H      9.8922015254   16.1733530537   -9.0391165650\n"
                            "H     11.4830844974   16.2742556871   -7.0270210879\n"
                            "H     12.3191478540   15.6022613957   -8.4131349788\n");
  std::stringstream struct2("60\n\n"
                            "C     13.9627782190   -2.7031343336   -3.6819011564\n"
                            "C     13.5930016397   -1.5051120752   -4.5473481240\n"
                            "C     12.1203553376   -1.5170526177   -4.9735633209\n"
                            "C     11.7522058777   -0.3111131073   -5.8434769880\n"
                            "C     10.2766646720   -0.3249915104   -6.2687436312\n"
                            "C      9.9270727450    0.8338372350   -7.1643163632\n"
                            "C      8.8286817381    1.5992996799   -7.0998744342\n"
                            "C      7.6975666859    1.4697765028   -6.1172451834\n"
                            "C      7.2505345618    2.8300203030   -5.5306042418\n"
                            "H      8.1265798237    3.3230256796   -5.0453355178\n"
                            "O      6.2605890576    2.5555801849   -4.5311580934\n"
                            "C      6.6739257527    3.7637488671   -6.5754379007\n"
                            "C      7.2744154009    4.8994449122   -6.9841715785\n"
                            "C      6.7076960993    5.7841460161   -7.9946459613\n"
                            "C      7.2446053295    6.9384006613   -8.4406320569\n"
                            "C      8.5307475987    7.5726597092   -7.9803315562\n"
                            "C      9.4361906889    7.9320845014   -9.1347584855\n"
                            "C     10.1019123462    9.0824519100   -9.3037025779\n"
                            "C     10.0950757770   10.2771824082   -8.3865723389\n"
                            "C      9.9898387783   11.6007019998   -9.1568693369\n"
                            "C     10.0313383814   12.8192526075   -8.2299713740\n"
                            "C      9.8705857622   14.1491816431   -8.9898172110\n"
                            "O      9.5270323929   14.1836264072  -10.1535412728\n"
                            "C     10.1926971190   15.4455090433   -8.2238889574\n"
                            "C      9.5834896165   15.5184390642   -6.8281362167\n"
                            "H     15.0118281830   -2.6667080208   -3.3958996684\n"
                            "H     13.7920326299   -3.6370543540   -4.2136870113\n"
                            "H     13.3671813777   -2.7262971302   -2.7714868998\n"
                            "H     14.2295746118   -1.4876682917   -5.4384302950\n"
                            "H     13.8060919032   -0.5800753383   -4.0010886107\n"
                            "H     11.4829871622   -1.5346302583   -4.0829105960\n"
                            "H     11.9057735816   -2.4420941397   -5.5198142778\n"
                            "H     12.3885841917   -0.2931297907   -6.7353081021\n"
                            "H     11.9668381362    0.6132918548   -5.2970588864\n"
                            "H      9.6367808922   -0.3372550776   -5.3810277894\n"
                            "H     10.0724105279   -1.2645324865   -6.8021621377\n"
                            "H     10.6521685363    1.0488375921   -7.9503565867\n"
                            "H      8.7169780845    2.3952649123   -7.8372918734\n"
                            "H      7.9673037449    0.8091336052   -5.2891139479\n"
                            "H      6.8326985789    1.0071274325   -6.6113209414\n"
                            "H      5.8642493788    3.3882160596   -4.2335907555\n"
                            "H      5.7137677377    3.4810394303   -7.0082849527\n"
                            "H      8.2341321952    5.1749239191   -6.5453182487\n"
                            "H      5.7565265454    5.4646886891   -8.4212120472\n"
                            "H      6.7005316013    7.4906973975   -9.2065942851\n"
                            "H      9.0645185185    6.8976353561   -7.2991136726\n"
                            "H      8.2967389794    8.4721802683   -7.3946997467\n"
                            "H      9.5578576232    7.1555771445   -9.8895495768\n"
                            "H     10.7315888106    9.1815335578  -10.1885341906\n"
                            "H     11.0236941626   10.2833510181   -7.7973818057\n"
                            "H      9.2784137121   10.2065164688   -7.6606845369\n"
                            "H      9.0622146250   11.6237480941   -9.7344005919\n"
                            "H     10.8051811282   11.6707846842   -9.8830329536\n"
                            "H     10.9661853116   12.8413118279   -7.6581602004\n"
                            "H      9.2295206422   12.7516500574   -7.4845836898\n"
                            "H      9.9606741440   14.7240663222   -6.1872935045\n"
                            "H      9.8599600940   16.2955803764   -8.8237312715\n"
                            "H     11.2835662654   15.5262545328   -8.1467937612\n"
                            "H      8.4999337007   15.4321002109   -6.8686890139\n"
                            "H      9.8250214716   16.4684296669   -6.3566030367\n"

  );

  PositionCollection pos1 = XyzStreamHandler::read(struct1).getPositions();
  PositionCollection pos2 = XyzStreamHandler::read(struct2).getPositions();
  indices = Distances::getListOfDivergingAtomsRobust(pos1, pos2, 0.05, 1e-2, 20);
  auto indices2 = Distances::getListOfDivergingAtomsRobust(pos1, pos2, 0.3, 1e-2, 20);
  for (int index : indices)
    ASSERT_TRUE(std::find(indices.begin(), indices.end(), index) != indices.end());
}

TEST(GeometryTest, GetsNearestNeighborsCorrectlyMethane) {
  std::stringstream methane("5\n\n"
                            "C      0.0000000000    0.0000000000    0.0000000000\n"
                            "H      0.6287000000    0.6287000000    0.6287000000\n"
                            "H     -0.6287000000   -0.6287000000    0.6287000000\n"
                            "H     -0.6287000000    0.6287000000   -0.6287000000\n"
                            "H      0.6287000000   -0.6287000000   -0.6287000000\n");
  auto structure = XyzStreamHandler::read(methane);
  BondOrderCollection bo = Geometry::Distances::nearestNeighborsBondOrders(structure.getPositions());
  EXPECT_EQ(bo.getMatrix().diagonal(), Eigen::VectorXd::Zero(bo.getSystemSize()));
  BondOrderCollection expected = BondOrderCollection(5);
  expected.setOrder(0, 1, 1.0);
  expected.setOrder(0, 2, 1.0);
  expected.setOrder(0, 3, 1.0);
  expected.setOrder(0, 4, 1.0);
  EXPECT_TRUE(bo.getMatrix().isApprox(expected.getMatrix()));
}

TEST(GeometryTest, GetsNearestNeighborsCorrectlyCH3Cl) {
  std::stringstream ch3cl("5\n\n"
                          "C        -0.08733982625534    0.08733981632076    0.08733982402444\n"
                          "H         0.55787462647471    0.70460025621186    0.70460025917197\n"
                          "H        -0.70460027220175   -0.55787461760986    0.70460027126431\n"
                          "H        -0.70460026285272    0.70460025895517   -0.55787462442182\n"
                          "Cl        0.93866573483510   -0.93866571387792   -0.93866573003890\n");
  auto structure = XyzStreamHandler::read(ch3cl);
  BondOrderCollection bo = Geometry::Distances::nearestNeighborsBondOrders(structure.getPositions());
  EXPECT_EQ(bo.getMatrix().diagonal(), Eigen::VectorXd::Zero(bo.getSystemSize()));
  BondOrderCollection expected = BondOrderCollection(5);
  expected.setOrder(0, 1, 1.0);
  expected.setOrder(0, 2, 1.0);
  expected.setOrder(0, 3, 1.0);
  expected.setOrder(0, 4, 1.0);
  EXPECT_TRUE(bo.getMatrix().isApprox(expected.getMatrix()));
}

TEST(GeometryTest, GetsNumberOfNearestNeighborsCorrectlyMethane) {
  std::stringstream methane("5\n\n"
                            "C      0.0000000000    0.0000000000    0.0000000000\n"
                            "H      0.6287000000    0.6287000000    0.6287000000\n"
                            "H     -0.6287000000   -0.6287000000    0.6287000000\n"
                            "H     -0.6287000000    0.6287000000   -0.6287000000\n"
                            "H      0.6287000000   -0.6287000000   -0.6287000000\n");
  auto structure = XyzStreamHandler::read(methane);
  auto nNeighbors = Geometry::Distances::countAllNearestNeighbors(structure.getPositions());
  EXPECT_EQ(nNeighbors.at(0), 4);
  for (unsigned i = 1; i < nNeighbors.size(); ++i) {
    EXPECT_EQ(nNeighbors.at(i), 1);
  }
}

TEST(GeometryTest, CanCalculateRotTransFreeTranformMatrix) {
  std::stringstream methane("5\n\n"
                            "C      0.0000000000    0.0000000000    0.0000000000\n"
                            "H      0.6287000000    0.6287000000    0.6287000000\n"
                            "H     -0.6287000000   -0.6287000000    0.6287000000\n"
                            "H     -0.6287000000    0.6287000000   -0.6287000000\n"
                            "H      0.6287000000   -0.6287000000   -0.6287000000\n");
  auto structure = XyzStreamHandler::read(methane);

  GradientCollection arbitraryGradient = GradientCollection::Zero(5, 3);
  arbitraryGradient.row(1) = Position{-1, -1, -1};

  Eigen::MatrixXd transformMatrix = Geometry::Transformations::calculateRotTransFreeTransformMatrix(
      structure.getPositions(), structure.getElements(), arbitraryGradient, true);

  Eigen::VectorXd overlaps =
      transformMatrix.transpose() * Eigen::Map<const Eigen::VectorXd>(arbitraryGradient.data(), arbitraryGradient.size());

  EXPECT_THAT(overlaps.norm(), DoubleNear(0.0, 1e-6));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
