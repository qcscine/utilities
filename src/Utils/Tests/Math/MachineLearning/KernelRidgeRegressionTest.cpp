/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Math/MachineLearning/CrossValidation.h>
#include <Utils/Math/MachineLearning/Regression/KernelRidgeRegression.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Utils {
using namespace MachineLearning;
using namespace testing;
namespace Tests {

/**
 * @class AKernelRidgeRegressionTest @file KernelRidgeRegressionTest.cpp
 * @test
 *
 * The KRR is tested against the implementation of scikit-learn:
 * https://scikit-learn.org/stable/modules/generated/sklearn.kernel_ridge.KernelRidge.html (visited Jan 14, 2020)
 *
 */
class AKernelRidgeRegressionTest : public Test {
 public:
  Eigen::MatrixXd linearDataFeatures;
  Eigen::VectorXd linearDataTargets;
  Eigen::MatrixXd nonlinearDataFeatures;
  Eigen::MatrixXd nonlinearDataTargets;

 private:
  void SetUp() final {
    // Set-up for the linear model
    // Underlying function: f(x,y,z) = 7.5x + 3.2y + 1.8z
    linearDataFeatures.resize(30, 3);
    for (int i = 0; i < 10; ++i) {
      linearDataFeatures(i, 0) = i;
      linearDataFeatures(i, 1) = i;
      linearDataFeatures(i, 2) = i;
    }
    for (int i = 10; i < 20; ++i) {
      linearDataFeatures(i, 0) = i;
      linearDataFeatures(i, 1) = i - 0.1;
      linearDataFeatures(i, 2) = i + 0.1;
    }
    for (int i = 20; i < 30; ++i) {
      linearDataFeatures(i, 0) = i + 0.1;
      linearDataFeatures(i, 1) = i + 0.2;
      linearDataFeatures(i, 2) = i - 0.2;
    }
    // Filling the targets with the correct values
    linearDataTargets.resize(30);
    for (int j = 0; j < 30; ++j) {
      linearDataTargets(j) = 7.5 * linearDataFeatures(j, 0) + 3.2 * linearDataFeatures(j, 1) + 1.8 * linearDataFeatures(j, 2);
    }

    // Set-up for the non-linear model
    // Underlying functions for the two targets
    // g(x) = x - 2 * x^3 * cos(x)
    // h(x) = x^2
    nonlinearDataFeatures.resize(10, 1);
    for (int i = 0; i < 10; ++i) {
      nonlinearDataFeatures(i, 0) = 0.5 * i * M_PI / 9;
    }
    nonlinearDataTargets.resize(10, 2);
    for (int j = 0; j < 10; ++j) {
      auto x = nonlinearDataFeatures(j, 0);
      nonlinearDataTargets(j, 0) = x - 2 * pow(x, 3) * cos(x);
      nonlinearDataTargets(j, 1) = x * x;
    }
  }
};

TEST_F(AKernelRidgeRegressionTest, LinearRegressionWorksCorrectly) {
  KernelRidgeRegression krr;
  krr.setRegularizationFactor(1e-4); // Small but non-zero
  krr.trainModel(linearDataFeatures, linearDataTargets);

  Eigen::VectorXd pointToPredict(3);
  pointToPredict(0) = 7.3;
  pointToPredict(1) = 7.1;
  pointToPredict(2) = 6.8;

  double p = krr.predict(pointToPredict)(0);
  EXPECT_THAT(p, DoubleNear(89.71, 1e-1));
}

TEST_F(AKernelRidgeRegressionTest, CrossValidationWorksCorrectly) {
  KernelRidgeRegression krr;
  krr.setRegularizationFactor(1.0);

  CrossValidation cv(krr, 10);
  cv.setRandomSeed(42);
  auto result = cv.evaluateRegressionModel(linearDataFeatures, linearDataTargets);
  EXPECT_THAT(result.first, DoubleNear(0.15, 1e-2));
  EXPECT_THAT(result.second, DoubleNear(0.05, 1e-2));

  // With a smaller regularization factor, the results should get very small
  krr.setRegularizationFactor(1e-6);
  result = cv.evaluateRegressionModel(linearDataFeatures, linearDataTargets);
  EXPECT_TRUE(result.first < 5e-4 && result.second < 5e-4);

  // Leave-one-out cross validation
  krr.setRegularizationFactor(1.0);
  cv.setNumberOfSubsets(30);
  result = cv.evaluateRegressionModel(linearDataFeatures, linearDataTargets);
  EXPECT_THAT(result.first, DoubleNear(0.138, 1e-2));
  EXPECT_THAT(result.second, DoubleNear(0.074, 1e-2));
}

TEST_F(AKernelRidgeRegressionTest, NonlinearRegressionWorksCorrectly) {
  KernelRidgeRegression krr;
  krr.setRegularizationFactor(1.0);

  // First try a linear model
  krr.trainModel(nonlinearDataFeatures, nonlinearDataTargets);

  Eigen::VectorXd pointToPredict(1);
  pointToPredict(0) = 1.0;

  auto predictions = krr.predict(pointToPredict);
  double p1 = predictions(0);
  double p2 = predictions(1);
  // wrong solution as linear model is not sufficient, tested again scikit-learn
  EXPECT_THAT(p1, DoubleNear(0.35, 1e-2));
  // Also insufficient model for the second target
  EXPECT_TRUE(p2 > 1.1);

  // Try the Gaussian kernel with bad sigma first
  krr.setRegularizationFactor(1e-6);
  krr.setKernel(Kernels::gaussianKernel, {10.0});
  krr.trainModel(nonlinearDataFeatures, nonlinearDataTargets);

  predictions = krr.predict(pointToPredict);
  // wrong solution as because of bad sigma
  EXPECT_THAT(predictions(0), DoubleNear(0.15, 1e-2));
  // For the second target the sigma is good enough, but not perfect
  EXPECT_THAT(predictions(1), DoubleNear(1.0, 1e-2));
  EXPECT_TRUE(predictions(1) > 1.005);

  // Try the Gaussian kernel with good sigma and even better regularization factor
  krr.setRegularizationFactor(1e-8);
  krr.setKernel(Kernels::gaussianKernel, {1.0});
  krr.trainModel(nonlinearDataFeatures, nonlinearDataTargets);

  predictions = krr.predict(pointToPredict);
  EXPECT_THAT(predictions(0), DoubleNear(-0.0806046, 1e-3));
  EXPECT_THAT(predictions(1), DoubleNear(1.0, 1e-4));

  // Small sigma will lead to overfitting
  krr.setRegularizationFactor(1e-8);
  krr.setKernel(Kernels::gaussianKernel, {0.01});
  krr.trainModel(nonlinearDataFeatures, nonlinearDataTargets);

  predictions = krr.predict(pointToPredict);
  EXPECT_THAT(predictions(0), DoubleNear(0.0, 1e-4));
  EXPECT_THAT(predictions(1), DoubleNear(0.0, 1e-4));

  // Large regularization factor will also yield bad results for the first target prediction
  krr.setRegularizationFactor(10.0);
  krr.setKernel(Kernels::gaussianKernel, {1.0});
  krr.trainModel(nonlinearDataFeatures, nonlinearDataTargets);

  predictions = krr.predict(pointToPredict);
  EXPECT_TRUE(predictions(0) > 0.0);

  // Test that the other kernels also work
  // Laplacian kernel:
  krr.setRegularizationFactor(1e-8);
  krr.setKernel(Kernels::laplacianKernel, {1.0});
  krr.trainModel(nonlinearDataFeatures, nonlinearDataTargets);
  predictions = krr.predict(pointToPredict);
  // Good results, but not as good as the Gaussian kernel
  EXPECT_THAT(predictions(0), DoubleNear(-0.0806046, 2e-2));
  EXPECT_THAT(predictions(1), DoubleNear(1.0, 5e-3));

  // Polynomial kernel of third degree:
  krr.setRegularizationFactor(1e-8);
  krr.setKernel(Kernels::polynomialKernel, {3, 1.0, 1.0});
  krr.trainModel(nonlinearDataFeatures, nonlinearDataTargets);
  predictions = krr.predict(pointToPredict);
  EXPECT_THAT(predictions(0), DoubleNear(-0.0806046, 5e-2));
  EXPECT_THAT(predictions(1), DoubleNear(1.0, 5e-3));
  // difference is larger than e.g., applying the Gaussian kernel
  EXPECT_TRUE(std::abs(predictions(0) + 0.0806046) > 1e-3);

  // 5-fold cross validation
  CrossValidation cv(krr, 5);
  cv.setRandomSeed(42);
  auto performanceThirdDegreePolynomial = cv.evaluateRegressionModel(nonlinearDataFeatures, nonlinearDataTargets);

  // Polynomial kernel of fifth degree (better results)
  krr.setRegularizationFactor(1e-8);
  krr.setKernel(Kernels::polynomialKernel, {5, 1.0, 1.0});
  krr.trainModel(nonlinearDataFeatures, nonlinearDataTargets);
  predictions = krr.predict(pointToPredict);
  EXPECT_THAT(predictions(0), DoubleNear(-0.0806046, 1e-3));
  EXPECT_THAT(predictions(1), DoubleNear(1.0, 1e-4));

  // 5-fold cross validation
  auto performanceFifthDegreePolynomial = cv.evaluateRegressionModel(nonlinearDataFeatures, nonlinearDataTargets);

  EXPECT_TRUE(performanceFifthDegreePolynomial.first < performanceThirdDegreePolynomial.first);
  EXPECT_TRUE(performanceFifthDegreePolynomial.second < performanceThirdDegreePolynomial.second);
  EXPECT_TRUE(performanceFifthDegreePolynomial.first < 1e-2);
  EXPECT_TRUE(performanceThirdDegreePolynomial.first < 2e-1);
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
