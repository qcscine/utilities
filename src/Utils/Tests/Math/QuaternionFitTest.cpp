/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Math/QuaternionFit.h>
#include <gmock/gmock.h>
#include <Eigen/Geometry>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class QuaternionFitTest QuaternionFitTest.cpp
 * @brief Comprises tests for the class Scine::Utils::QuaternionFit.
 * @test
 */
class QuaternionFitTest : public Test {
 public:
  const int numberParticles = 22;
  Eigen::MatrixXd randomPos1, randomPos2;

 protected:
  void SetUp() override {
    randomPos1 = Eigen::MatrixXd::Random(numberParticles, 3) * 22;
    randomPos2 = Eigen::MatrixXd::Random(numberParticles, 3) * 10.6;
  }
};

TEST_F(QuaternionFitTest, ReturnsZeroTranslationForIdenticalMatrix) {
  QuaternionFit fit(randomPos1, randomPos1);
  Eigen::Vector3d translationVector = fit.getTransVector();
  ASSERT_TRUE(translationVector.isZero());
}

TEST_F(QuaternionFitTest, TwoDimers) {
  PositionCollection positions1(2, 3);
  PositionCollection positions2(2, 3);
  positions1 << -2, 0, 0, 2, 0, 0;
  positions2 << -1, 0, 0, 1, 0, 0;
  QuaternionFit fit(positions1, positions2);
  ASSERT_NEAR(fit.getRMSD(), 1.0, 1E-6);
}

TEST_F(QuaternionFitTest, Chirality) {
  std::stringstream stream("5\n\n"
                           "C     -0.0220968652    0.0032150545    0.0165197402\n"
                           "H     -0.6690087809    0.8893598645   -0.1009085033\n"
                           "H     -0.3777879439   -0.8577518853   -0.5882960277\n"
                           "H      0.0964209243   -0.3151252987    1.0637808670\n"
                           "H      0.9724726657    0.2803022650   -0.3910960761");

  auto atoms1 = XyzStreamHandler::read(stream);
  auto atoms2 = atoms1;

  // Swap positions
  atoms2.setPosition(1, atoms1.getPosition(2));
  atoms2.setPosition(2, atoms1.getPosition(1));

  // left hand vs right hand do not allow inversion (default), expect no inversion and bad fit
  QuaternionFit fit1(atoms1.getPositions(), atoms2.getPositions());
  ASSERT_NEAR(fit1.getRMSD(), 2.098277, 1E-6);

  Eigen::MatrixXd refMat = atoms1.getPositions();
  Eigen::MatrixXd fitMat = atoms2.getPositions();
  Eigen::VectorXd weights = Eigen::VectorXd::Ones(atoms1.size());
  bool improperRotationsAllowed = true;

  // left hand vs right hand allow inversion, expect inversion
  QuaternionFit fit2(refMat, fitMat, weights, improperRotationsAllowed);
  ASSERT_NEAR(fit2.getRMSD(), 0.046643, 1E-6);

  // left hand vs left hand allow inverison, expect none
  refMat = atoms1.getPositions();
  fitMat = atoms2.getPositions();
  QuaternionFit fit3(refMat, fitMat, weights, improperRotationsAllowed);
  ASSERT_NEAR(fit3.getRMSD(), 0.046643, 1E-6);
}

TEST_F(QuaternionFitTest, Weights) {
  std::stringstream stream("3\n\n"
                           "O    -0.00392892   0.38587356  0.00000000\n"
                           "H    -0.76036657  -0.19836384  0.00000000\n"
                           "H     0.75900101  -0.18964299  0.00000000");

  auto atoms1 = XyzStreamHandler::read(stream);
  auto atoms2 = atoms1;

  // Move one position
  auto pos = atoms1.getPosition(1);
  atoms2.setPosition(2, pos + Position{0, 0, 1});

  Eigen::MatrixXd refMat = atoms1.getPositions();
  Eigen::MatrixXd fitMat = atoms2.getPositions();
  Eigen::VectorXd fitWeights = Eigen::VectorXd::Ones(atoms1.size());
  fitWeights[2] = 0;
  bool improperRotationsAllowed = false;
  QuaternionFit fit(refMat, fitMat, fitWeights, improperRotationsAllowed);

  // Standard RMSD
  ASSERT_NEAR(fit.getRMSD(), 1.755372, 1E-6);

  // Weighted RMSD
  ASSERT_NEAR(fit.getWeightedRMSD(fitWeights), 0, 1E-6);
}

TEST_F(QuaternionFitTest, ReturnsZeroRotationForIdenticalMatrix) {
  QuaternionFit fit(randomPos1, randomPos1);
  Eigen::Matrix3d rotationMatrix = fit.getRotationMatrix();
  ASSERT_TRUE(rotationMatrix.isApprox(Eigen::Matrix3d::Identity()));
}

TEST_F(QuaternionFitTest, IsCorrectForTranslatedPositions) {
  Eigen::Vector3d translation(-9.654, 5.85588, 0.11);

  Eigen::MatrixXd translatedPositions = randomPos1;
  translatedPositions.rowwise() += translation.transpose();
  QuaternionFit fit(randomPos1, translatedPositions);

  Eigen::Vector3d translationVector = fit.getTransVector();
  Eigen::Matrix3d rotationMatrix = fit.getRotationMatrix();

  ASSERT_TRUE(translationVector.isApprox(translation));
  ASSERT_TRUE(rotationMatrix.isApprox(Eigen::Matrix3d::Identity()));
}

TEST_F(QuaternionFitTest, IsCorrectForRotatedPositions) {
  // First, necessary to center the reference positions:
  Eigen::Vector3d mean = randomPos1.colwise().mean();
  randomPos1.rowwise() -= mean.transpose();

  Eigen::AngleAxisd rotation(4.001, Eigen::Vector3d::Random().normalized());

  Eigen::MatrixXd rotatedPositions(numberParticles, 3);
  for (int i = 0; i < numberParticles; ++i) {
    rotatedPositions.row(i) = rotation * randomPos1.row(i).transpose();
  }
  QuaternionFit fit(randomPos1, rotatedPositions);

  Eigen::Vector3d translationVector = fit.getTransVector();
  Eigen::Matrix3d rotationMatrix = fit.getRotationMatrix();

  ASSERT_TRUE(translationVector.isZero());
  ASSERT_TRUE(rotationMatrix.isApprox(rotation.toRotationMatrix()));
}

} // namespace Tests
} /* namespace Utils */
} /* namespace Scine */
