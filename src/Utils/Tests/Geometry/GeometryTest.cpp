/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry.h"
#include "Utils/Math/QuaternionFit.h"
#include "Utils/Typenames.h"
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

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
