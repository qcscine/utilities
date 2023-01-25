/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Constants.h>
#include <Utils/DataStructures/PeriodicBoundaries.h>
#include <Utils/Geometry/Utilities/Distances.h>
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Utils {
namespace Tests {

class APeriodicBoundariesTest : public Test {
 public:
  PositionCollection randomPositions1;
  PositionCollection randomPositions2;
  PositionCollection refCartPositions1;
  PositionCollection refRelPositions1;
  PositionCollection refCartPositions2;
  PositionCollection refRelPositions2;
  double a1;
  double b1;
  double c1;
  double a2;
  double b2;
  double c2;
  double alpha1;
  double beta1;
  double gamma1;
  double alpha2;
  double beta2;
  double gamma2;
  Eigen::Vector3d lengths1;
  Eigen::Vector3d lengths2;
  Eigen::Vector3d angles1;
  Eigen::Vector3d angles2;

 private:
  // void TestBody() override {
  void SetUp() override {
    randomPositions1 = Eigen::MatrixX3d::Random(50, 3);
    randomPositions2 = Eigen::MatrixX3d::Random(50, 3);
    // all following data is in Angstrom and Degrees
    // cell 1
    a1 = 6.52372159;
    b1 = 6.52372159;
    c1 = 6.52372159;
    alpha1 = 90.00000000;
    beta1 = 90.00000000;
    gamma1 = 90.00000000;
    lengths1 << a1, b1, c1;
    angles1 << alpha1, beta1, gamma1;
    // cell 2
    a2 = 3.297078325;
    b2 = 3.29707832413;
    c2 = 5.254212582;
    alpha2 = 90.0;
    beta2 = 90.0;
    gamma2 = 120.0;
    lengths2 << a2, b2, c2;
    angles2 << alpha2, beta2, gamma2;
    // clang-format off
      // taken from https://github.com/materialsproject/pymatgen/blob/master/test_files/CuCl.cif
      refRelPositions1 = Eigen::MatrixX3d::Zero(16, 3);
      refRelPositions1 << 0.157170, 0.342830, 0.657170,
                          0.342830, 0.657170, 0.157170,
                          0.657170, 0.157170, 0.342830,
                          0.842830, 0.842830, 0.842830,
                          0.842830, 0.657170, 0.342830,
                          0.657170, 0.342830, 0.842830,
                          0.342830, 0.842830, 0.657170,
                          0.157170, 0.157170, 0.157170,
                          0.619186, 0.880814, 0.119186,
                          0.880814, 0.119186, 0.619186,
                          0.119186, 0.619186, 0.880814,
                          0.380814, 0.380814, 0.380814,
                          0.380814, 0.119186, 0.880814,
                          0.119186, 0.880814, 0.380814,
                          0.880814, 0.380814, 0.119186,
                          0.619186, 0.619186, 0.619186;
      refCartPositions1 = Eigen::MatrixX3d::Zero(16, 3);
      refCartPositions1 <<  1.02533, 2.23653, 4.28719,
                            2.23653, 4.28719, 1.02533,
                            4.28719, 1.02533, 2.23653,
                            5.49839, 5.49839, 5.49839,
                            5.49839, 4.28719, 2.23653,
                            4.28719, 2.23653, 5.49839,
                            2.23653, 5.49839, 4.28719,
                            1.02533, 1.02533, 1.02533,
                            4.03940, 5.74619, 0.77754,
                            5.74619, 0.77754, 4.03940,
                            0.77754, 4.03940, 5.74619,
                            2.48432, 2.48432, 2.48432,
                            2.48432, 0.77754, 5.74619,
                            0.77754, 5.74619, 2.48432,
                            5.74619, 2.48432, 0.77754,
                            4.03940, 4.03940, 4.03940;
      // taken from https://github.com/materialsproject/pymatgen/blob/master/test_files/CoO19128.cif
      refRelPositions2 = Eigen::MatrixX3d::Zero(4, 3);
      refRelPositions2 <<  0.666666127557, 0.333331747190,  0.878675905847,
                           0.333333498186, 0.666666996573,  0.378675199937,
                           0.666666127557, 0.333331747190,  0.496323656362,
                           0.333333498186, 0.666666996573,  0.996324362271;
      refCartPositions2 = Eigen::MatrixX3d::Zero(4, 3);
      refCartPositions2 <<  1.64854,       0.95178,       4.61675,
                           -0.00000,       1.90357,       1.98964,
                            1.64854,       0.95178,       2.60779,
                           -0.00000,       1.90357,       5.23490;
    // clang-format on
  }
};

TEST_F(APeriodicBoundariesTest, ClassIsInitializedCorrectly) {
  PeriodicBoundaries pbc0 = PeriodicBoundaries();
  pbc0 *= 6.52372159 * Constants::bohr_per_angstrom;
  PeriodicBoundaries pbc1 = PeriodicBoundaries(lengths1, angles1, false);
  EXPECT_TRUE(pbc0 == pbc1);
  EXPECT_DOUBLE_EQ(pbc1.getA().norm(), a1 * Constants::bohr_per_angstrom);
  EXPECT_DOUBLE_EQ(pbc1.getB().norm(), b1 * Constants::bohr_per_angstrom);
  EXPECT_DOUBLE_EQ(pbc1.getC().norm(), c1 * Constants::bohr_per_angstrom);
  EXPECT_DOUBLE_EQ(pbc1.getAlpha(), alpha1);
  EXPECT_DOUBLE_EQ(pbc1.getBeta(), beta1);
  EXPECT_DOUBLE_EQ(pbc1.getGamma(), gamma1);
  std::string expectedString = std::to_string(a1 * Constants::bohr_per_angstrom) + "," +
                               std::to_string(b1 * Constants::bohr_per_angstrom) + "," +
                               std::to_string(c1 * Constants::bohr_per_angstrom) + "," + std::to_string(alpha1) + "," +
                               std::to_string(beta1) + "," + std::to_string(gamma1) + ",xyz";
  EXPECT_EQ(pbc1.getPeriodicBoundariesString(), expectedString);

  PeriodicBoundaries pbc2 = PeriodicBoundaries(lengths2, angles2, false);
  EXPECT_DOUBLE_EQ(pbc2.getA().norm(), a2 * Constants::bohr_per_angstrom);
  EXPECT_DOUBLE_EQ(pbc2.getB().norm(), b2 * Constants::bohr_per_angstrom);
  EXPECT_DOUBLE_EQ(pbc2.getC().norm(), c2 * Constants::bohr_per_angstrom);
  EXPECT_DOUBLE_EQ(pbc2.getAlpha(), alpha2);
  EXPECT_DOUBLE_EQ(pbc2.getBeta(), beta2);
  EXPECT_DOUBLE_EQ(pbc2.getGamma(), gamma2);

  Eigen::Vector3d lengths3 = lengths1 * Constants::bohr_per_angstrom;
  PeriodicBoundaries pbc3 = PeriodicBoundaries(lengths3, angles1);
  EXPECT_DOUBLE_EQ(pbc3.getA().norm(), a1 * Constants::bohr_per_angstrom);
  EXPECT_DOUBLE_EQ(pbc3.getB().norm(), b1 * Constants::bohr_per_angstrom);
  EXPECT_DOUBLE_EQ(pbc3.getC().norm(), c1 * Constants::bohr_per_angstrom);
  EXPECT_DOUBLE_EQ(pbc3.getAlpha(), alpha1);
  EXPECT_DOUBLE_EQ(pbc3.getBeta(), beta1);
  EXPECT_DOUBLE_EQ(pbc3.getGamma(), gamma1);
  EXPECT_EQ(pbc1.getCellMatrix(), pbc3.getCellMatrix());

  Eigen::Vector3d angles4 = angles2;
  angles4[0] *= Constants::pi / 180.0;
  angles4[1] *= Constants::pi / 180.0;
  angles4[2] *= Constants::pi / 180.0;
  PeriodicBoundaries pbc4 = PeriodicBoundaries(lengths2, angles4, false, false);
  EXPECT_DOUBLE_EQ(pbc4.getA().norm(), a2 * Constants::bohr_per_angstrom);
  EXPECT_DOUBLE_EQ(pbc4.getB().norm(), b2 * Constants::bohr_per_angstrom);
  EXPECT_DOUBLE_EQ(pbc4.getC().norm(), c2 * Constants::bohr_per_angstrom);
  EXPECT_DOUBLE_EQ(pbc4.getAlpha(), alpha2);
  EXPECT_DOUBLE_EQ(pbc4.getBeta(), beta2);
  EXPECT_DOUBLE_EQ(pbc4.getGamma(), gamma2);
  EXPECT_EQ(pbc4.getCellMatrix(), pbc2.getCellMatrix());
  EXPECT_EQ(pbc4.getA(), pbc2.getA());
  EXPECT_EQ(pbc4.getB(), pbc2.getB());
  EXPECT_EQ(pbc4.getC(), pbc2.getC());

  Eigen::Matrix3d matrix = pbc4.getCellMatrix();
  PeriodicBoundaries pbc5 = PeriodicBoundaries(matrix);
  EXPECT_EQ(pbc5, pbc4);

  std::string input = std::to_string(a2) + ";" + std::to_string(b2) + "; " + std::to_string(c2) + ";" +
                      std::to_string(alpha2) + "; " + std::to_string(beta2) + ";" + std::to_string(gamma2);
  PeriodicBoundaries pbc6 = PeriodicBoundaries(input, ";", false);
  EXPECT_TRUE(pbc2.getCellMatrix().isApprox(pbc6.getCellMatrix(), 1e-5));
  input += ";xz";
  PeriodicBoundaries pbc7 = PeriodicBoundaries(input, ";", false);
  EXPECT_TRUE(pbc7 != pbc6);
  ASSERT_THAT(pbc7.getPeriodicityString(), "xz");
  std::array<bool, 3> expected = {true, false, true};
  ASSERT_THAT(pbc7.getPeriodicity(), expected);

  // reassignments
  pbc2.setCellMatrix(pbc1.getCellMatrix());
  ASSERT_TRUE(pbc1 == pbc2);
  ASSERT_DOUBLE_EQ(pbc1.getA().norm(), pbc2.getA().norm());
  ASSERT_DOUBLE_EQ(pbc1.getAlpha(), pbc2.getAlpha());
  pbc1.setPeriodicity("xy");
  pbc4 = pbc1;
  ASSERT_TRUE(pbc1 == pbc4);
  ASSERT_DOUBLE_EQ(pbc1.getA().norm(), pbc4.getA().norm());
  ASSERT_DOUBLE_EQ(pbc1.getAlpha(), pbc4.getAlpha());
  PeriodicBoundaries pbc8 = PeriodicBoundaries(Eigen::Matrix3d::Identity(), pbc1.getPeriodicityString());
  pbc1 = Eigen::Matrix3d::Identity();
  ASSERT_TRUE(pbc1 == pbc8);
}

TEST_F(APeriodicBoundariesTest, OrthoRhombicWorks) {
  PeriodicBoundaries pbc1 = PeriodicBoundaries(lengths1, angles1, false);
  PeriodicBoundaries pbc2 = PeriodicBoundaries(lengths2, angles2, false);
  ASSERT_TRUE(pbc1.isOrthoRhombic());
  ASSERT_FALSE(pbc2.isOrthoRhombic());
}

TEST_F(APeriodicBoundariesTest, CanonicalizationWorks) {
  Eigen::Matrix3d randomMatrix = Eigen::Matrix3d::Random(3, 3);
  for (int i = 0; i < 3; ++i) {
    randomMatrix(i, i) = std::fabs(randomMatrix(i, i));
  }
  auto pbc = PeriodicBoundaries(randomMatrix);
  auto canonic = PeriodicBoundaries(pbc.getPeriodicBoundariesString());
  auto rot = pbc.getCanonicalizationRotationMatrix();
  Eigen::Matrix3d manualMatrix = pbc.getCellMatrix() * rot;
  ASSERT_TRUE(std::fabs(canonic.getCellMatrix()(0, 2) < 1e-12));
  ASSERT_TRUE(std::fabs(canonic.getCellMatrix()(1, 2) < 1e-12));
  ASSERT_TRUE(std::fabs(manualMatrix(0, 2) < 1e-12));
  ASSERT_TRUE(std::fabs(manualMatrix(1, 2) < 1e-12));
  ASSERT_TRUE(canonic.getCellMatrix().isApprox(manualMatrix, 1e-6));
  pbc.canonicalize();
  ASSERT_TRUE(canonic.getCellMatrix().isApprox(pbc.getCellMatrix(), 1e-6));
  ASSERT_TRUE(std::fabs(pbc.getCellMatrix()(0, 2) < 1e-12));
  ASSERT_TRUE(std::fabs(pbc.getCellMatrix()(1, 2) < 1e-12));
  ASSERT_TRUE(pbc.getCanonicalizationRotationMatrix().isApprox(Eigen::Matrix3d::Identity()));
}

TEST_F(APeriodicBoundariesTest, LengthsAndAngles) {
  Eigen::Vector3d lengths = lengths2 * Constants::bohr_per_angstrom;
  auto pbc = PeriodicBoundaries(lengths, angles2);
  ASSERT_THAT(pbc.getLengths().size(), 3);
  ASSERT_THAT(pbc.getAngles().size(), 3);
  ASSERT_DOUBLE_EQ(pbc.getLengths()[0], lengths[0]);
  ASSERT_DOUBLE_EQ(pbc.getLengths()[1], lengths[1]);
  ASSERT_DOUBLE_EQ(pbc.getLengths()[2], lengths[2]);
  ASSERT_DOUBLE_EQ(pbc.getAngles()[0], angles2[0]);
  ASSERT_DOUBLE_EQ(pbc.getAngles()[1], angles2[1]);
  ASSERT_DOUBLE_EQ(pbc.getAngles()[2], angles2[2]);
}

TEST_F(APeriodicBoundariesTest, InverseIsConsistent) {
  PeriodicBoundaries pbc = PeriodicBoundaries(lengths1, angles1, false);
  ASSERT_TRUE(pbc.getInverseCellMatrix().isApprox(pbc.getCellMatrix().inverse()));
}

TEST_F(APeriodicBoundariesTest, MultiplicationsWork) {
  auto pbc1 = PeriodicBoundaries(10.0);
  auto pbc2 = PeriodicBoundaries(20.0);
  auto pbc3 = pbc1 * 2.0;
  auto scaling = Eigen::Vector3d::Constant(2.0);
  auto pbc4 = pbc1 * scaling;
  pbc1 *= 2.0;
  pbc1 *= 0.5;
  pbc1 *= scaling;
  EXPECT_EQ(pbc1, pbc2);
  EXPECT_EQ(pbc1, pbc3);
  EXPECT_EQ(pbc1, pbc4);
  EXPECT_EQ(pbc2, pbc3);
  EXPECT_EQ(pbc2, pbc4);
  Eigen::Vector3d scaling2;
  scaling2 << 10.0, 2.5, 4.0;
  pbc1 *= scaling2;
  EXPECT_DOUBLE_EQ(pbc1.getA().norm(), 200.0);
  EXPECT_DOUBLE_EQ(pbc1.getB().norm(), 50.0);
  EXPECT_DOUBLE_EQ(pbc1.getC().norm(), 80.0);
  EXPECT_DOUBLE_EQ(pbc1.getAlpha(), 90.0);
}

TEST_F(APeriodicBoundariesTest, IsSelfConsistent) {
  PeriodicBoundaries pbc1 = PeriodicBoundaries(lengths1, angles1, false);
  PositionCollection cartPositions1 = pbc1.transform(refCartPositions1, false);
  pbc1.transformInPlace(cartPositions1);
  EXPECT_TRUE(cartPositions1.isApprox(refCartPositions1, 1e-12));

  PositionCollection pos = refRelPositions1.row(2);
  PositionCollection relPosition = pbc1.transform(pos);
  pbc1.transformInPlace(relPosition, false);
  EXPECT_TRUE(relPosition.isApprox(pos, 1e-12));

  PeriodicBoundaries pbc2 = PeriodicBoundaries(lengths2, angles2, false);
  PositionCollection relPositions2 = pbc2.transform(refRelPositions2);
  pbc2.transformInPlace(relPositions2, false);
  EXPECT_TRUE(relPositions2.isApprox(refRelPositions2, 1e-12));

  pos = refCartPositions2.row(3);
  pbc2.transformInPlace(pos, false);
  PositionCollection cartPosition = pbc2.transform(pos);
  EXPECT_TRUE(cartPosition.isApprox(refCartPositions2.row(3), 1e-12));
}

TEST_F(APeriodicBoundariesTest, RelativeToCartesianIsCorrect) {
  PeriodicBoundaries pbc1 = PeriodicBoundaries(lengths1, angles1, false);
  PositionCollection cart1 = pbc1.transform(refRelPositions1);
  PositionCollection refCartBohr = refCartPositions1 * Constants::bohr_per_angstrom;
  EXPECT_TRUE(cart1.isApprox(refCartBohr, 1e-6));

  PeriodicBoundaries pbc2 = PeriodicBoundaries(lengths2, angles2, false);
  PositionCollection cart2 = pbc2.transform(refRelPositions2);
  refCartBohr = refCartPositions2 * Constants::bohr_per_angstrom;
  EXPECT_TRUE(cart2.isApprox(refCartBohr, 1e-6));
}

TEST_F(APeriodicBoundariesTest, CartesianToRelativeIsCorrect) {
  PeriodicBoundaries pbc1 = PeriodicBoundaries(lengths1, angles1, false);
  PositionCollection refCartBohr = refCartPositions1 * Constants::bohr_per_angstrom;
  PositionCollection rel1 = pbc1.transform(refCartBohr, false);
  EXPECT_TRUE(rel1.isApprox(refRelPositions1, 1e-6));

  PeriodicBoundaries pbc2 = PeriodicBoundaries(lengths2, angles2, false);
  refCartBohr = refCartPositions2 * Constants::bohr_per_angstrom;
  PositionCollection rel2 = pbc2.transform(refCartBohr, false);
  EXPECT_TRUE(rel2.isApprox(refRelPositions2, 1e-6));
}

TEST_F(APeriodicBoundariesTest, CanTranslateIntoCell) {
  PeriodicBoundaries pbc1 = PeriodicBoundaries(lengths1, angles1, false);
  PeriodicBoundaries pbc2 = PeriodicBoundaries(lengths2, angles2, false);
  Position farthestPos = Eigen::RowVector3d::Ones();
  double biggestDistInCell1 = (pbc1.transform(farthestPos)).norm();
  double biggestDistInCell2 = (pbc2.transform(farthestPos)).norm();
  EXPECT_LT(pbc1.translatePositionsIntoCell(pbc1.transform(farthestPos)).norm(), biggestDistInCell1);
  EXPECT_LT(pbc2.translatePositionsIntoCell(pbc2.transform(farthestPos)).norm(), biggestDistInCell2);

  randomPositions1 *= 9.11;
  PositionCollection pos = pbc1.translatePositionsIntoCell(randomPositions1);
  for (int i = 0; i < randomPositions1.rows(); ++i) {
    Position relPos = pos.row(i);
    EXPECT_LE(relPos.norm(), biggestDistInCell1);
    pbc1.transformInPlace(relPos, false);
    EXPECT_LE(relPos[0], 1.0);
    EXPECT_GE(relPos[0], 0.0);
    EXPECT_LE(relPos[1], 1.0);
    EXPECT_GE(relPos[1], 0.0);
    EXPECT_LE(relPos[2], 1.0);
    EXPECT_GE(relPos[2], 0.0);
  }
  pbc2.translatePositionsIntoCellInPlace(randomPositions1);
  for (int i = 0; i < randomPositions1.rows(); ++i) {
    Position relPos = randomPositions1.row(i);
    EXPECT_LE(relPos.norm(), biggestDistInCell2);
    pbc2.transformInPlace(relPos, false);
    EXPECT_LE(relPos[0], 1.0);
    EXPECT_GE(relPos[0], 0.0);
    EXPECT_LE(relPos[1], 1.0);
    EXPECT_GE(relPos[1], 0.0);
    EXPECT_LE(relPos[2], 1.0);
    EXPECT_GE(relPos[2], 0.0);
  }
}

TEST_F(APeriodicBoundariesTest, WithinCheckWorks) {
  PeriodicBoundaries pbcHuge = PeriodicBoundaries(Eigen::Matrix3d::Identity() * 100.0);
  pbcHuge.translatePositionsIntoCellInPlace(randomPositions1);
  ASSERT_TRUE(pbcHuge.isWithinCell(randomPositions1));
  PeriodicBoundaries pbcSmall = PeriodicBoundaries(Eigen::Matrix3d::Identity() * 2.0);
  PositionCollection positions = Eigen::MatrixX3d::Zero(2, 3);
  positions << 0.1, 1.1, 1.7, 1.0, 1e-12, 2.3;
  ASSERT_FALSE(pbcSmall.isWithinCell(positions));
  ASSERT_TRUE(pbcSmall.isWithinCell(static_cast<Position>(positions.row(0))));
  ASSERT_FALSE(pbcSmall.isWithinCell(static_cast<Position>(positions.row(1))));
  pbcSmall.translatePositionsIntoCellInPlace(positions);
  ASSERT_TRUE(pbcSmall.isWithinCell(positions));
}

TEST_F(APeriodicBoundariesTest, CalculatesDistanceCorrectly) {
  PeriodicBoundaries pbc1 = PeriodicBoundaries(lengths1, angles1, false);
  PeriodicBoundaries pbc2 = PeriodicBoundaries(lengths2, angles2, false);
  Position farthestPos = Eigen::RowVector3d::Ones();
  Position farthestCartesian1 = (pbc1.transform(farthestPos));
  Position farthestCartesian2 = (pbc2.transform(farthestPos));
  Position origin = Eigen::RowVector3d::Zero();
  EXPECT_NEAR(Geometry::Distances::distance(origin, farthestCartesian1, pbc1), 0.0, 1e-12);
  EXPECT_NEAR(Geometry::Distances::distance(origin, farthestCartesian2, pbc2), 0.0, 1e-12);

  Position middlePos;
  middlePos << 0.5, 0.5, 0.5;
  Position middleCartesian1 = pbc1.transform(middlePos);
  Position middleCartesian2 = pbc2.transform(middlePos);
  double biggestMinimumImageDist1 = std::sqrt(pbc1.getBiggestPossibleDistanceSquared());
  double biggestMinimumImageDist2 = std::sqrt(pbc2.getBiggestPossibleDistanceSquared());
  EXPECT_NEAR(Geometry::Distances::distance(origin, middleCartesian1, pbc1), biggestMinimumImageDist1, 1e-12);
  EXPECT_NEAR(Geometry::Distances::distance(origin, middleCartesian2, pbc2), biggestMinimumImageDist2, 1e-12);

  int nAtoms = randomPositions1.rows();
  EXPECT_LE(Geometry::Distances::distance(randomPositions1, randomPositions2, pbc1), biggestMinimumImageDist1 * nAtoms);
  EXPECT_LE(Geometry::Distances::distance(randomPositions1, randomPositions2, pbc2), biggestMinimumImageDist2 * nAtoms);
  for (int i = 0; i < nAtoms; ++i) {
    Position pos1 = randomPositions1.row(i);
    Position pos2 = randomPositions2.row(i);
    EXPECT_LE(Geometry::Distances::distance(pos1, pos2, pbc1), biggestMinimumImageDist1);
    EXPECT_LE(Geometry::Distances::distance(pos1, pos2, pbc2), biggestMinimumImageDist2);
  }

  Position pos;
  pos << 0.0, 0.0, 0.2;
  Position cartPos1 = pbc1.transform(pos);
  Position cartPos2 = pbc2.transform(pos);
  double dist1 = pbc1.getC().norm() * 0.2;
  double dist2 = pbc2.getC().norm() * 0.2;
  EXPECT_NEAR(Geometry::Distances::distance(origin, cartPos1, pbc1), dist1, 1e-12);
  EXPECT_NEAR(Geometry::Distances::distance(origin, cartPos2, pbc2), dist2, 1e-12);

  pos(2) = 0.8;
  cartPos1 = pbc1.transform(pos);
  cartPos2 = pbc2.transform(pos);
  EXPECT_NEAR(Geometry::Distances::distance(origin, cartPos1, pbc1), dist1, 1e-12);
  EXPECT_NEAR(Geometry::Distances::distance(origin, cartPos2, pbc2), dist2, 1e-12);
}

TEST_F(APeriodicBoundariesTest, CalculatesMinimumDisplacementCorrectly) {
  PeriodicBoundaries pbc1 = PeriodicBoundaries(lengths1, angles1, false);
  Position pos1;
  pos1 << 0.0, 0.2, 0.0;
  Position pos2;
  pos2 << 0.0, 0.0, 0.2;
  pbc1.translatePositionsIntoCellInPlace(pos1);
  pbc1.translatePositionsIntoCellInPlace(pos1);
  Displacement disp1 = pbc1.bruteForceMinimumImageDisplacementVector(pos1, pos2);
  Displacement disp2 = pbc1.bruteForceMinimumImageDisplacementVector(pos2, pos1);
  ASSERT_TRUE(disp1.isApprox(pos2 - pos1));
  ASSERT_TRUE(disp2.isApprox(pos1 - pos2));
  // pos3 in cell below should be closer to pos1 than same cell
  Position pos3;
  pos3 << 0.1, 0.0, c1 * Constants::bohr_per_angstrom - 0.3;
  pbc1.translatePositionsIntoCellInPlace(pos3);
  Displacement disp3Exp;
  disp3Exp << 0.1, -0.2, -0.3;
  Displacement disp3 = pbc1.bruteForceMinimumImageDisplacementVector(pos1, pos3);
  ASSERT_TRUE(disp3.isApprox(disp3Exp));
  Displacement disp4 = pbc1.bruteForceMinimumImageDisplacementVector(pos3, pos1);
  ASSERT_TRUE(disp4.isApprox(-disp3Exp));
}

TEST_F(APeriodicBoundariesTest, GetsNearestNeighborsCorrectly) {
  PositionCollection cart1 = refCartPositions1 * Constants::bohr_per_angstrom;
  PositionCollection cart2 = refCartPositions2 * Constants::bohr_per_angstrom;
  PositionCollection cart3 = Eigen::MatrixX3d::Zero(4, 3);
  // clang-format off
  // FM-3M space group unit cell taken from https://materialsproject.org/materials/mp-30/
  cart3 << 0.00000000, 0.00000000, 0.00000000,
           0.00000000, 0.50000000, 0.50000000,
           0.50000000, 0.00000000, 0.50000000,
           0.50000000, 0.50000000, 0.00000000;
  // clang-format on
  Eigen::Vector3d lengths;
  Eigen::Vector3d angles;
  lengths << 2.561, 2.561, 2.561;
  angles << 90.0, 90.0, 90.0;

  PeriodicBoundaries pbc1 = PeriodicBoundaries(lengths1, angles1, false);
  PeriodicBoundaries pbc2 = PeriodicBoundaries(lengths2, angles2, false);
  PeriodicBoundaries pbc3 = PeriodicBoundaries(lengths, angles, false);
  pbc3.transformInPlace(cart3);

  BondOrderCollection bo1 = Geometry::Distances::nearestNeighborsBondOrders(cart1, pbc1);
  BondOrderCollection bo2 = Geometry::Distances::nearestNeighborsBondOrders(cart2, pbc2);
  BondOrderCollection bo3 = Geometry::Distances::nearestNeighborsBondOrders(cart3, pbc3);

  EXPECT_EQ(bo1.getMatrix().diagonal(), Eigen::VectorXd::Zero(bo1.getSystemSize()));
  EXPECT_EQ(bo2.getMatrix().diagonal(), Eigen::VectorXd::Zero(bo2.getSystemSize()));
  Eigen::MatrixXd allBonded = Eigen::MatrixXd::Ones(4, 4);
  allBonded.diagonal() = Eigen::VectorXd::Zero(4);
  EXPECT_TRUE(bo3.getMatrix().isApprox(allBonded));

  BondOrderCollection expected1 = BondOrderCollection(16);
  // clang-format off
  expected1.setOrder(9,  0, 1.0);
  expected1.setOrder(10, 0, 1.0);
  expected1.setOrder(11, 0, 1.0);
  expected1.setOrder(8,  1, 1.0);
  expected1.setOrder(10, 1, 1.0);
  expected1.setOrder(11, 1, 1.0);
  expected1.setOrder(8,  2, 1.0);
  expected1.setOrder(9,  2, 1.0);
  expected1.setOrder(11, 2, 1.0);
  expected1.setOrder(8,  3, 1.0);
  expected1.setOrder(9,  3, 1.0);
  expected1.setOrder(10, 3, 1.0);
  expected1.setOrder(13, 4, 1.0);
  expected1.setOrder(14, 4, 1.0);
  expected1.setOrder(15, 4, 1.0);
  expected1.setOrder(12, 5, 1.0);
  expected1.setOrder(14, 5, 1.0);
  expected1.setOrder(15, 5, 1.0);
  expected1.setOrder(12, 6, 1.0);
  expected1.setOrder(13, 6, 1.0);
  expected1.setOrder(15, 6, 1.0);
  expected1.setOrder(14, 7, 1.0);
  expected1.setOrder(13, 7, 1.0);
  expected1.setOrder(12, 7, 1.0);
  // clang-format on
  EXPECT_TRUE(bo1.getMatrix().isApprox(expected1.getMatrix()));

  Eigen::MatrixXd expected2 = Eigen::MatrixXd::Zero(4, 4);
  // clang-format off
  expected2 << 0.0, 0.0, 1.0, 1.0,
          0.0, 0.0, 1.0, 1.0,
          1.0, 1.0, 0.0, 0.0,
          1.0, 1.0, 0.0, 0.0;
  // clang-format on
  EXPECT_TRUE(Eigen::MatrixXd(bo2.getMatrix()).isApprox(expected2));
}
TEST_F(APeriodicBoundariesTest, GetsNumberOfNearestNeighborsCorrectly) {
  PositionCollection cart1 = refCartPositions1 * Constants::bohr_per_angstrom;
  PositionCollection cart2 = refCartPositions2 * Constants::bohr_per_angstrom;
  PositionCollection cart3 = Eigen::MatrixX3d::Zero(4, 3);
  // clang-format off
  // FM-3M space group unit cell taken from https://materialsproject.org/materials/mp-30/
  cart3 << 0.00000000, 0.00000000, 0.00000000,
           0.00000000, 0.50000000, 0.50000000,
           0.50000000, 0.00000000, 0.50000000,
           0.50000000, 0.50000000, 0.00000000;
  // clang-format on
  Eigen::Vector3d lengths;
  Eigen::Vector3d angles;
  lengths << 2.561, 2.561, 2.561;
  angles << 90.0, 90.0, 90.0;

  PeriodicBoundaries pbc1 = PeriodicBoundaries(lengths1, angles1, false);
  PeriodicBoundaries pbc2 = PeriodicBoundaries(lengths2, angles2, false);
  PeriodicBoundaries pbc3 = PeriodicBoundaries(lengths, angles, false);
  pbc3.transformInPlace(cart3);

  auto nNeighbors1 = Geometry::Distances::countAllNearestNeighbors(cart1, pbc1);
  auto nNeighbors2 = Geometry::Distances::countAllNearestNeighbors(cart2, pbc2);
  auto nNeighbors3 = Geometry::Distances::countAllNearestNeighbors(cart3, pbc3);

  for (const auto& n : nNeighbors1) {
    EXPECT_EQ(n, 3);
  }

  for (const auto& n : nNeighbors2) {
    EXPECT_EQ(n, 4);
  }

  for (const auto& n : nNeighbors3) {
    EXPECT_EQ(n, 12);
  }
}

TEST_F(APeriodicBoundariesTest, GetsNumberOfNearestNeighborsOfOneCorrectly) {
  Eigen::Vector3d lengths;
  Eigen::Vector3d angles;
  lengths << 2.0, 2.0, 2.0;
  angles << 90.0, 90.0, 90.0;
  PeriodicBoundaries pbc = PeriodicBoundaries(lengths, angles);

  PositionCollection coord = Eigen::MatrixX3d::Zero(2, 3);
  coord << 0.5, 0.5, 0.5, 0.0, 0.0, 0.0;
  pbc.transformInPlace(coord);
  EXPECT_EQ(Geometry::Distances::countNearestNeighbors(coord, 0, pbc), 8);
  EXPECT_EQ(Geometry::Distances::countNearestNeighbors(coord, 1, pbc), 8);
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
