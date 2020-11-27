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
  translatedTestPC = translatePositions(testPC, shiftVector);

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

  auto matrix = positionVectorToMatrix(positionVector);

  for (int i = 0; i < numberParticles; ++i) {
    for (int j = 0; j < 3; ++j) {
      ASSERT_THAT(matrix(i, j), DoubleEq(positionVector(3 * i + j)));
    }
  }
}

TEST(GeometryTest, CanTransformAPositionMatrixToAVector) {
  int numberParticles = 11;

  Eigen::MatrixXd positionMatrix = Eigen::MatrixXd::Random(numberParticles, 3);

  auto vector = positionMatrixToVector(positionMatrix);

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

  PositionCollection rotatedPC = rotatePositions(testPC, startOrientation, endOrientation, rotOrigin);
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
  translatePositionsInPlace(waterPC0, axisOrg1);

  // Rotation around z-axis of 120 degrees
  double rotAngle1 = 120 * Constants::pi / 180;
  Displacement axis1(0, 0, 1);

  PositionCollection waterPC1 = rotatePositions(waterPC0, axis1, rotAngle1, axisOrg1);
  Position v1 = waterPC0.row(1) - waterPC0.row(0);
  Position v2 = waterPC1.row(1) - waterPC1.row(0);
  double dotPr1 = v1.dot(v2);
  double measuredAngle1 = acos(dotPr1 / (v1.norm() * v2.norm())) * 180 / Constants::pi;

  // Rotation around x-axis of 90 degrees
  double rotAngle2 = 90 * Constants::pi / 180;
  Displacement axis2(1, 0, 0);

  PositionCollection waterPC2 = rotatePositions(waterPC0, axis2, rotAngle2, axisOrg1);
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
  randomDisplacementInPlace(pc1, maxDisplacement);
  for (int i = 0; i < pc1.rows(); ++i) {
    for (int j = 0; j < pc1.cols(); ++j) {
      ASSERT_LE(std::fabs(pc1(i, j) - pc2(i, j)), maxDisplacement);
    }
  }
}

TEST(GeometryTest, CanGenerateRandomTrajectoryFromAtomCollection) {
  PositionCollection pc(2, 3);
  pc << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
  auto ec = ElementTypeCollection{ElementType::H, ElementType::H};
  AtomCollection atoms(ec, pc);
  int numFrames = 10;
  double maxDisplacement = 0.1;
  MolecularTrajectory trajectory = randomDisplacementTrajectory(atoms, numFrames, maxDisplacement);
  EXPECT_THAT(trajectory.molecularSize(), pc.rows());
  EXPECT_THAT(trajectory.size(), numFrames);
  for (int i = 0; i < trajectory.size(); ++i) {
    for (int j = 0; j < pc.rows(); ++j) {
      for (int k = 0; k < pc.cols(); ++k) {
        EXPECT_LE(std::fabs(trajectory[i](j, k) - pc(j, k)), maxDisplacement);
      }
    }
  }
}

TEST(GeometryTest, CanAlignTwoPositionCollections) {
  int numberParticles = 11;
  PositionCollection referencePositions = Eigen::MatrixXd::Random(numberParticles, 3);
  PositionCollection positions = Eigen::MatrixXd::Random(numberParticles, 3);

  alignPositions(referencePositions, positions);

  QuaternionFit fit(referencePositions, positions);

  Eigen::Vector3d translationVector = fit.getTransVector();
  Eigen::Matrix3d rotation = fit.getRotationMatrix();

  Eigen::AngleAxisd angleRotation{rotation};

  ASSERT_THAT(angleRotation.angle(), DoubleNear(0.0, 1e-6));
  ASSERT_TRUE(translationVector.isZero(1e-6));
}

TEST(GeometryTest, GetsCorrectMassVector) {
  auto ec = ElementTypeCollection{ElementType::H, ElementType::C, ElementType::Fe};

  auto masses = getMasses(ec);
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

  auto index = getIndexOfAtomInStructure(structure, atom1);
  ASSERT_THAT(index, Eq(2));

  Position p2(-0.6287, -0.7287, 0.6287);
  p2 *= Constants::bohr_per_angstrom;
  Atom atom2(ElementType::H, p2);
  Atom atom3(ElementType::Cl, p1);

  EXPECT_THROW(getIndexOfAtomInStructure(structure, atom2), std::runtime_error);
  EXPECT_THROW(getIndexOfAtomInStructure(structure, atom3), std::runtime_error);

  Position p3(-0.6287, -0.63, 0.6287);
  p3 *= Constants::bohr_per_angstrom;
  Atom atom4(ElementType::H, p3);

  EXPECT_THROW(getIndexOfAtomInStructure(structure, atom3), std::runtime_error);
  ASSERT_THAT(getIndexOfAtomInStructure(structure, atom4, 0.1), Eq(2));
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

  ASSERT_THAT(getIndexOfClosestAtom(positions, p1), Eq(2));
  ASSERT_THAT(getIndexOfClosestAtom(positions, p1, -0.001), Eq(2));
  ASSERT_THAT(getIndexOfClosestAtom(positions, p1, 1e-4), Eq(0));

  Position p2(-1.0, 1.0, -1.0);
  p2 *= Constants::bohr_per_angstrom;

  ASSERT_THAT(getIndexOfClosestAtom(positions, p2), Eq(3));
  ASSERT_THAT(getIndexOfClosestAtom(positions, p2, -0.001), Eq(3));
  ASSERT_THAT(getIndexOfClosestAtom(positions, p2, 1e-4), Eq(3));
}

TEST(GeometryTest, Centroid) {
  PositionCollection diametric(2, 3);
  diametric << 1.0, 2.0, 3.0, -1.0, 2.0, -3.0;
  ASSERT_TRUE(getAveragePosition(diametric).isApprox(Position(0.0, 2.0, 0.0), 1e-10));
}

TEST(GeometryTest, CanRecoverListOfDivergingAtomPositions) {
  PositionCollection reference(5, 3);
  reference << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0, 4.0;
  PositionCollection testPositions = reference;
  EXPECT_EQ(getListOfDivergingAtoms(reference, testPositions, 0.11).size(), 0);

  DisplacementCollection deviation = DisplacementCollection::Zero(5, 3);
  deviation(3, 2) -= 0.5;
  deviation(4, 2) += 0.5;
  testPositions += deviation;

  auto indices = getListOfDivergingAtoms(reference, testPositions, 0.11);
  EXPECT_EQ(indices.size(), 2);
  EXPECT_EQ(indices[0], 3);
  EXPECT_EQ(indices[1], 4);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
