/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/StructuralCompletion.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {
/**
 * @class StructuralCompletionTest StructuralCompletionTest.cpp
 * @brief Comprises tests for the class Scine::Utils::AtomCollection.
 * @test
 */
class StructuralCompletionTest : public Test {
 public:
  AtomCollection methane;
  static void assertAngleIsCorrect(const Eigen::RowVector3d v1, const Eigen::RowVector3d v2, double correctAngle);
  double tetrahedronAngle;
  double triangularAngle;

 protected:
  void SetUp() override {
    tetrahedronAngle = 109.5;
    triangularAngle = 120.0;
    std::stringstream methaneCoordinates("5\n\n"
                                         "H     -0.6287000000   -0.6287000000    0.6287000000\n"
                                         "H     -0.6287000000    0.6287000000   -0.6287000000\n"
                                         "C      0.0000000000    0.0000000000    0.0000000000\n"
                                         "H      0.6287000000   -0.6287000000   -0.6287000000\n"
                                         "H      0.6287000000    0.6287000000    0.6287000000\n");
    methane = XyzStreamHandler::read(methaneCoordinates);
  }
};

TEST_F(StructuralCompletionTest, ThreeTetrahedralCornersAreBuiltFromOneOther) {
  std::vector<Eigen::Vector3d> newPositions(3);
  // get distance vector from center of tetraedron to corner
  Eigen::Vector3d v1 = (methane.at(1).getPosition() - methane.at(2).getPosition()).normalized();
  StructuralCompletion::generate3TetrahedronCornersFrom1Other(v1, newPositions.at(0), newPositions.at(1), newPositions.at(2));

  // new distance vectors should be normalized
  for (auto& position : newPositions) {
    ASSERT_THAT(position.norm(), DoubleNear(1.00, 1e-3));
  }

  // angle between distance vectors should correspond to tetrahedral angle
  for (unsigned int i = 0; i < newPositions.size() - 1; i++) {
    assertAngleIsCorrect(newPositions.at(i), newPositions.at(i + 1), tetrahedronAngle);
  }
}

TEST_F(StructuralCompletionTest, TwoTetrahedralCornersAreBuiltFromTwoOthers) {
  std::vector<Eigen::Vector3d> newPositions(2);
  Eigen::Vector3d v1 = (methane.at(1).getPosition() - methane.at(2).getPosition()).normalized();
  Eigen::Vector3d v2 = (methane.at(0).getPosition() - methane.at(2).getPosition()).normalized();

  StructuralCompletion::generate3TetrahedronCornersFrom1Other(v1, v2, newPositions.at(0), newPositions.at(1));

  for (auto& position : newPositions) {
    ASSERT_THAT(position.norm(), DoubleNear(1.00, 1e-3));
  }

  assertAngleIsCorrect(v1, newPositions.at(0), tetrahedronAngle);
  assertAngleIsCorrect(newPositions.at(0), newPositions.at(1), tetrahedronAngle);
}

TEST_F(StructuralCompletionTest, OneTetrahedralCornerIsBuiltFromThreeOthers) {
  Eigen::Vector3d newPos;
  Eigen::Vector3d v1 = (methane.at(1).getPosition() - methane.at(2).getPosition()).normalized();
  Eigen::Vector3d v2 = (methane.at(0).getPosition() - methane.at(2).getPosition()).normalized();
  Eigen::Vector3d v3 = (methane.at(3).getPosition() - methane.at(2).getPosition()).normalized();

  StructuralCompletion::generate3TetrahedronCornersFrom1Other(v1, v2, v3, newPos);

  ASSERT_THAT(newPos.norm(), DoubleNear(1.00, 1e-3));

  assertAngleIsCorrect(newPos, v1, tetrahedronAngle);
  assertAngleIsCorrect(newPos, v2, tetrahedronAngle);
  assertAngleIsCorrect(newPos, v3, tetrahedronAngle);
}

TEST_F(StructuralCompletionTest, TwoTriangleCornersAreBuiltFromOneOther) {
  std::vector<Eigen::Vector3d> newPositions(2);
  Eigen::Vector3d v1 = {1.0, 0.0, 0.0};
  Eigen::Vector3d mid = {0.0, 0.0, 0.0};
  Eigen::Vector3d dist1 = (v1 - mid).normalized();

  StructuralCompletion::generate2TriangleCornersFrom1Other(dist1, newPositions.at(0), newPositions.at(1));

  ASSERT_THAT(newPositions.at(0).norm(), DoubleNear(1.00, 1e-3));
  ASSERT_THAT(newPositions.at(1).norm(), DoubleNear(1.00, 1e-3));

  assertAngleIsCorrect(newPositions.at(0), dist1, triangularAngle);
  assertAngleIsCorrect(newPositions.at(1), dist1, triangularAngle);
}

TEST_F(StructuralCompletionTest, OneTriangleCornerIsBuiltFromTwoOthers) {
  Eigen::Vector3d newPos;
  Eigen::Vector3d v1 = {1.0, 0.0, 0.0};
  Eigen::Vector3d mid = {0.0, 0.0, 0.0};
  Eigen::Vector3d v2 = {-0.5, -0.866025, 0.0};

  Eigen::Vector3d dist1 = (v1 - mid).normalized();
  Eigen::Vector3d dist2 = (v2 - mid).normalized();

  StructuralCompletion::generate1TriangleCornerFrom2Others(dist1, dist2, newPos);

  ASSERT_THAT(newPos.norm(), DoubleNear(1.00, 1e-3));

  assertAngleIsCorrect(newPos, dist1, triangularAngle);
  assertAngleIsCorrect(newPos, dist2, triangularAngle);
}

void StructuralCompletionTest::assertAngleIsCorrect(const Eigen::RowVector3d v1, const Eigen::RowVector3d v2,
                                                    double correctAngle) {
  double angle = (Utils::Constants::degree_per_rad * acos(v1.dot(v2) / (v1.norm() * v2.norm())));
  ASSERT_THAT(angle, DoubleNear(correctAngle, 1e-1));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
