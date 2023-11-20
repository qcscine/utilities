/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Math/MachineLearning/Regression/GaussianProcessRegression.h>
#include <Utils/Math/MachineLearning/Regression/RegressionSettings.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Utils {
using namespace MachineLearning;
using namespace testing;
namespace Tests {

/**
 * @class AGaussianProcessRegressionTest @file AGaussianProcessRegressionTest.cpp
 * @test
 *
 */
class AGaussianProcessRegressionTest : public Test {
 public:
  Eigen::MatrixXd linearDataFeatures;
  Eigen::VectorXd linearDataTargets;

  Eigen::MatrixXd featureValues;
  Eigen::MatrixXd targetValues;

 private:
  void SetUp() final {
    featureValues.resize(7, 1);
    // X
    featureValues.col(0) << -4, -3, -2, -1, 1, 2, 3;
    // Y
    targetValues.resize(7, 1);
    targetValues = featureValues.array().sin();
  }
};

TEST_F(AGaussianProcessRegressionTest, FeaturesAndTargetsHaveCorrectDimension) {
  Eigen::MatrixXd featureValues(5, 1);
  featureValues.col(0) << -4, -3, -2, -1, 0;
  Eigen::MatrixXd targetValues(4, 1);
  targetValues << 0, 1, 2, 3;

  GaussianProcessRegression gpr;
  ASSERT_THROW(gpr.trainModel(featureValues, targetValues), std::runtime_error);
}

TEST_F(AGaussianProcessRegressionTest, SettingsCanBeModified) {
  GaussianProcessRegression gpr;
  gpr.settings().modifyInt(Utils::SettingsNames::Optimizations::MachineLearning::maxIterations, 10);
  gpr.settings().modifyInt(Utils::SettingsNames::Optimizations::MachineLearning::maxLinesearch, 100);
  gpr.settings().modifyDouble(Utils::SettingsNames::Optimizations::MachineLearning::convergenceTolerance, 1e-2);
  gpr.settings().modifyDouble(Utils::SettingsNames::Optimizations::MachineLearning::ftol, 1e-3);

  ASSERT_THAT(gpr.settings().getInt(Utils::SettingsNames::Optimizations::MachineLearning::maxIterations), Eq(10));
  ASSERT_THAT(gpr.settings().getInt(Utils::SettingsNames::Optimizations::MachineLearning::maxLinesearch), Eq(100));
  ASSERT_THAT(gpr.settings().getDouble(Utils::SettingsNames::Optimizations::MachineLearning::convergenceTolerance),
              DoubleNear(1e-2, 1e-2));
  ASSERT_THAT(gpr.settings().getDouble(Utils::SettingsNames::Optimizations::MachineLearning::ftol), DoubleNear(1e-3, 1e-3));
}

TEST_F(AGaussianProcessRegressionTest, SinusIsLearnedCorrectly) {
  GaussianProcessRegression gpr;
  // no noise
  gpr.setSigmaYSqHyperparametersGuess({0.0, false, {{0.0, 5}}});
  gpr.setSigmaFSqHyperparametersGuess({1.0, true, {{-1000, 1000}}});
  gpr.setThetaHyperparametersGuess({1.0, true, {{-1000, 1000}}});
  gpr.trainModel(featureValues, targetValues);

  Eigen::VectorXd testSet(34, 1);
  testSet << -3.80, -3.60, -3.40, -3.20, -3.00, -2.80, -2.60, -2.40, -2.20, -2.00, -1.80, -1.60, -1.40, -1.20, -1.00,
      -0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80;
  auto mean = gpr.predict(testSet);
  auto hyperparam = gpr.getOptimizedHyperparameters();
  // calculate the correct solution
  auto correctSolutions = testSet.array().sin();
  for (int i = 0; i < mean.size(); i++) {
    ASSERT_THAT(mean[i], DoubleNear(correctSolutions[i], 1e-2));
  };

  // with noise
  gpr.setSigmaYSqHyperparametersGuess({1.0, true, {{0.1, 10.0}}});
  gpr.setSigmaFSqHyperparametersGuess({1.0, true, {{0.0, std::numeric_limits<double>::infinity()}}});
  gpr.setThetaHyperparametersGuess({1.0, true, {{-100, 100.0}}});
  gpr.trainModel(featureValues, targetValues);
  Eigen::VectorXd testSet2(8, 1);
  // add a data point that is far out of the range of training points
  testSet2 << -3.80, -2.80, -1.80, -0.80, 0.80, 1.80, 2.80, 10.0;
  auto mean2 = gpr.predict(testSet2);
  auto covariance = gpr.getVarianceOfPrediction();
  ASSERT_TRUE(covariance.isApprox(covariance.transpose(), 1e-2));
  // check if variance is calculated correctly
  auto confidenceInterval = 1.96 * covariance.diagonal().array().sqrt();

  ASSERT_THAT(confidenceInterval[0], DoubleNear(0.50622676, 1e-4));
  ASSERT_THAT(confidenceInterval[1], DoubleNear(0.49089664, 1e-4));
  // high uncertainty for the last data point
  ASSERT_THAT(confidenceInterval[7], DoubleNear(2.1375408, 1e-4));
}

TEST_F(AGaussianProcessRegressionTest, LinearFunctionIsLearnedCorrectly) {
  // Set-up for the linear model
  // Underlying function: f(x,y,z) = 7.5x + 3.2y + 1.8z
  linearDataFeatures.resize(11, 1);
  for (int i = 0; i < 11; ++i) {
    linearDataFeatures(i, 0) = i;
  }
  // Filling the targets with the correct values
  linearDataTargets.resize(11);
  for (int j = 0; j < 11; ++j) {
    linearDataTargets(j) = 7.5 * linearDataFeatures(j, 0) + 3.2;
  }
  GaussianProcessRegression gpr;
  gpr.setSigmaYSqHyperparametersGuess({0.01, false, {{0.0, 0.0}}});
  gpr.setSigmaFSqHyperparametersGuess({10.0, true, {{-100, 100}}});
  gpr.setThetaHyperparametersGuess({10.0, true, {{-1000, 1000}}});

  gpr.trainModel(linearDataFeatures, linearDataTargets);
  Eigen::VectorXd testSet(12, 1);
  testSet << -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5;
  Eigen::VectorXd correctSolution(12);
  for (int i = 0; i < 12; i++) {
    testSet(i, 0) = 2.0 + 0.2 * i;
    correctSolution[i] = 7.5 * testSet(i, 0) + 3.2;
  }

  auto mean = gpr.predict(testSet);
  for (int i = 0; i < mean.size(); i++) {
    ASSERT_THAT(mean[i], DoubleNear(correctSolution[i], 1e-2));
  }
}

TEST_F(AGaussianProcessRegressionTest, HyperParameterAreCorrectlyOptimized) {
  GaussianProcessRegression gpr;
  // optimize only Theta with boundaries
  gpr.setSigmaYSqHyperparametersGuess({0.0, false, {{0.0, 0.0}}});
  gpr.setSigmaFSqHyperparametersGuess({1.0, false, {{1.0, 1.0}}});
  gpr.setThetaHyperparametersGuess({1.0, true, {{100, -100}}});
  ASSERT_THROW(gpr.trainModel(featureValues, targetValues), std::runtime_error);

  // optimize only Theta with boundaries
  gpr.setSigmaYSqHyperparametersGuess({0.0, false, {{0.0, 0.0}}});
  gpr.setSigmaFSqHyperparametersGuess({1.0, false, {{1.0, 1.0}}});
  gpr.setThetaHyperparametersGuess({1.0, true, {{-100, 100}}});

  gpr.trainModel(featureValues, targetValues);
  auto optParams1 = gpr.getOptimizedHyperparameters();
  ASSERT_THAT(optParams1[0], DoubleNear(1.60, 1e-2));
  ASSERT_THAT(optParams1[1], DoubleNear(1.0, 1e-2));
  ASSERT_THAT(optParams1[2], DoubleNear(0.0, 1e-2));

  // optimize only SigmaF
  gpr.setSigmaYSqHyperparametersGuess({0.1, false, {{0.1, 0.1}}});
  gpr.setSigmaFSqHyperparametersGuess({1.0, true, {{0, std::numeric_limits<double>::infinity()}}});
  gpr.setThetaHyperparametersGuess({1.0, false, {{1.0, 1.0}}});
  gpr.trainModel(featureValues, targetValues);
  auto optParams2 = gpr.getOptimizedHyperparameters();
  ASSERT_THAT(optParams2[0], DoubleNear(1.0, 1e-2));
  ASSERT_THAT(optParams2[1], DoubleNear(0.912, 1e-2));
  ASSERT_THAT(optParams2[2], DoubleNear(0.1, 1e-2));

  // optimize only SigmaY^2
  gpr.setSigmaYSqHyperparametersGuess({1.0, true, {{0.1, 10.0}}});
  gpr.setSigmaFSqHyperparametersGuess({1.0, false, {{1.0, 1.0}}});
  gpr.setThetaHyperparametersGuess({1.0, false, {{1.0, 1.0}}});
  gpr.trainModel(featureValues, targetValues);
  auto optParams3 = gpr.getOptimizedHyperparameters();
  ASSERT_THAT(optParams3[0], DoubleNear(1.0, 1e-2));
  ASSERT_THAT(optParams3[1], DoubleNear(1.0, 1e-2));
  ASSERT_THAT(optParams3[2], DoubleNear(0.1, 1e-2));

  // optimize all
  gpr.setSigmaYSqHyperparametersGuess({1.0, true, {{0.1, 10.0}}});
  gpr.setSigmaFSqHyperparametersGuess({1.0, true, {{0, std::numeric_limits<double>::infinity()}}});
  gpr.setThetaHyperparametersGuess({1.0, true, {{-100, 100}}});

  auto boundaries = gpr.getBoundaries();
  ASSERT_THAT(boundaries[0].first, Eq(0));
  ASSERT_THAT(boundaries[0].second.first, DoubleNear(-100.0, 1e-4));
  ASSERT_THAT(boundaries[0].second.second, DoubleNear(100.0, 1e-4));

  ASSERT_THAT(boundaries[1].first, Eq(1));
  ASSERT_THAT(boundaries[1].second.first, DoubleNear(0.0, 1e-4));
  ASSERT_THAT(boundaries[1].second.second, DoubleNear(std::numeric_limits<double>::infinity(), 1e-4));

  ASSERT_THAT(boundaries[2].first, Eq(2));
  ASSERT_THAT(boundaries[2].second.first, DoubleNear(0.1, 1e-4));
  ASSERT_THAT(boundaries[2].second.second, DoubleNear(10.0, 1e-4));

  gpr.trainModel(featureValues, targetValues);

  auto optParams4 = gpr.getOptimizedHyperparameters();
  ASSERT_THAT(optParams4[0], DoubleNear(1.409, 1e-2));
  ASSERT_THAT(optParams4[1], DoubleNear(1.091, 1e-2));
  ASSERT_THAT(optParams4[2], DoubleNear(0.10, 1e-2));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
