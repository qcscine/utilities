/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Optimizer/LeastSquares/LevenbergMarquardt.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

class UpdateFunctionMock : public UpdateFunctionManagerBase {
 public:
  UpdateFunctionMock(const Eigen::VectorXd& xValues, const Eigen::VectorXd& yValues)
    : xValues_(xValues), yValues_(yValues){};
  void updateErrors(const Eigen::VectorXd& parameters, Eigen::VectorXd& errors) override {
    int numDataPoints = getNumberOfDataPoints(parameters);
    errors.resize(numDataPoints);
    for (int i = 0; i < numDataPoints; ++i) {
      errors(i) = yValues_(i) - func(parameters, xValues_(i));
    }
  };

  int getNumberOfDataPoints(const Eigen::VectorXd& /*parameters*/) const override {
    return xValues_.size();
  }

 private:
  // f(x) = ax^3 + bx^2 + cx + d
  double func(const Eigen::VectorXd& parameters, double x) {
    return parameters(0) * std::pow(x, 3) + parameters(1) * std::pow(x, 2) + parameters(2) * x + parameters(3);
  }
  const Eigen::VectorXd& xValues_;
  const Eigen::VectorXd& yValues_;
};

class ALevenbergMarquardtTest : public Test {
 public:
  Eigen::VectorXd startingParameters;
  Eigen::VectorXd xValues;
  Eigen::VectorXd yValues;

 protected:
  void SetUp() override {
    startingParameters.resize(4);
    startingParameters(0) = 1.0;
    startingParameters(1) = 1.0;
    startingParameters(2) = 1.0;
    startingParameters(3) = 1.0;
    xValues.resize(6);
    yValues.resize(6);
    xValues(0) = -5.1;
    yValues(0) = 4.4;
    xValues(1) = -4.3;
    yValues(1) = -2.9;
    xValues(2) = -1.2;
    yValues(2) = 0.69;
    xValues(3) = 0.6;
    yValues(3) = 6.0;
    xValues(4) = 2.7;
    yValues(4) = -7.85;
    xValues(5) = 2.0;
    yValues(5) = 0.26;
  }
};

// The 6 data points were taken very close to, but not quite on, the actual curve
// f(x) = -0.4x^3 -1.6x^2 + 2.4x + 5.2
TEST_F(ALevenbergMarquardtTest, APolynomialOfThirdDegreeIsFittedReasonablyWell) {
  LevenbergMarquardt lm;
  UpdateFunctionMock update(xValues, yValues);
  lm.calculateCovarianceMatrix = true;
  lm.optimize(startingParameters, update);
  auto covarianceMatrix = lm.getCovarianceMatrix();

  ASSERT_THAT(startingParameters(0), DoubleNear(-0.4, 5e-2));
  ASSERT_THAT(startingParameters(1), DoubleNear(-1.6, 5e-2));
  ASSERT_THAT(startingParameters(2), DoubleNear(2.4, 5e-2));
  ASSERT_THAT(startingParameters(3), DoubleNear(5.2, 5e-2));
  for (int i = 0; i < startingParameters.size(); ++i) {
    for (int j = 0; j < startingParameters.size(); ++j) {
      ASSERT_TRUE(std::abs(covarianceMatrix(i, j)) > 1e-7);
    }
  }
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
