/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Optimizer/GradientBased/Bofill.h"
#include "Utils/Optimizer/GradientBased/LBFGS.h"
#include "Utils/Optimizer/GradientBased/SteepestDescent.h"
#include "Utils/Optimizer/HessianBased/EigenvectorFollowing.h"
#include "Utils/Optimizer/HessianBased/NewtonRaphson.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Settings.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

auto const gradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
  value = -cos(parameters[0]) - 0.5 * cos(parameters[1]);
  gradients[0] = sin(parameters[0]);
  gradients[1] = 0.5 * sin(parameters[1]);
};

auto const hessianTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                    Eigen::MatrixXd& hessian) {
  value = -cos(parameters[0]) - 0.5 * cos(parameters[1]);
  gradients[0] = sin(parameters[0]);
  gradients[1] = 0.5 * sin(parameters[1]);
  hessian(0, 0) = cos(parameters[0]);
  hessian(0, 1) = 0.0;
  hessian(1, 0) = 0.0;
  hessian(1, 1) = 0.5 * cos(parameters[1]);
};

/**
 * @class Scine::Utils::Tests::OptimizerTests
 * @brief Comprises tests for the class Scine::Utils::Optimizer derived classes.
 * @test
 */
TEST(OptimizerTests, SDTest) {
  SteepestDescent optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.25 * M_PI;
  positions[1] = 0.75 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 1000;
  check.deltaValue = 1e-12;
  check.stepMaxCoeff = 1.0e-7;
  check.stepRMS = 1.0e-9;
  check.gradMaxCoeff = 1.0e-7;
  check.gradRMS = 1.0e-9;
  check.requirement = 4;
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check);
  EXPECT_TRUE(nCycles < 1000);
  EXPECT_NEAR(positions[0], 0.0, 2.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 2.0e-8);
}
TEST(OptimizerTests, LBFGSTest) {
  LBFGS optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.25 * M_PI;
  positions[1] = 0.75 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 100;
  check.deltaValue = 1e-12;
  check.stepMaxCoeff = 1.0e-7;
  check.stepRMS = 1.0e-8;
  check.gradMaxCoeff = 1.0e-7;
  check.gradRMS = 1.0e-8;
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, NewtonRaphsonTest_Minimization) {
  NewtonRaphson optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.25 * M_PI;
  positions[1] = 0.49 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  check.stepMaxCoeff = 1.0e-7;
  check.stepRMS = 1.0e-8;
  check.gradMaxCoeff = 1.0e-7;
  check.gradRMS = 1.0e-8;
  auto nCycles = optimizer.optimize(positions, hessianTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, NewtonRaphsonTest_PartialMaximization) {
  NewtonRaphson optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.25 * M_PI;
  positions[1] = 0.51 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  check.stepMaxCoeff = 1.0e-7;
  check.stepRMS = 1.0e-8;
  check.gradMaxCoeff = 1.0e-7;
  check.gradRMS = 1.0e-8;
  auto nCycles = optimizer.optimize(positions, hessianTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], M_PI, 1.0e-8);
}
TEST(OptimizerTests, NewtonRaphsonTest_Maximization) {
  NewtonRaphson optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.51 * M_PI;
  positions[1] = 0.51 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, hessianTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], M_PI, 1.0e-8);
}
TEST(OptimizerTests, EigenvectorFollowingTest_1) {
  EigenvectorFollowing optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.49 * M_PI;
  positions[1] = 0.51 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, hessianTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], M_PI, 1.0e-8);
}
TEST(OptimizerTests, EigenvectorFollowingTest_2) {
  EigenvectorFollowing optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.51 * M_PI;
  positions[1] = 0.49 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, hessianTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, EigenvectorFollowingTest_3) {
  EigenvectorFollowing optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.51 * M_PI;
  positions[1] = 0.51 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, hessianTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, BofillTest_1) {
  Bofill optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.55 * M_PI;
  positions[1] = 0.55 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 100;
  check.deltaValue = 1e-12;
  check.stepMaxCoeff = 1.0e-7;
  check.stepRMS = 1.0e-8;
  check.gradMaxCoeff = 1.0e-7;
  check.gradRMS = 1.0e-8;
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
}
TEST(OptimizerTests, BofillTest_2) {
  Bofill optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.49 * M_PI;
  positions[1] = 0.51 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 100;
  check.deltaValue = 1e-12;
  check.stepMaxCoeff = 1.0e-7;
  check.stepRMS = 1.0e-8;
  check.gradMaxCoeff = 1.0e-7;
  check.gradRMS = 1.0e-8;
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
}
TEST(OptimizerTests, BofillTest_3) {
  Bofill optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.55 * M_PI;
  positions[1] = 0.45 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 100;
  check.deltaValue = 1e-12;
  check.stepMaxCoeff = 1.0e-7;
  check.stepRMS = 1.0e-8;
  check.gradMaxCoeff = 1.0e-7;
  check.gradRMS = 1.0e-8;
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
}
TEST(OptimizerTests, BofillTest_4) {
  Bofill optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.90 * M_PI;
  positions[1] = 0.10 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 100;
  check.deltaValue = 1e-12;
  check.stepMaxCoeff = 1.0e-7;
  check.stepRMS = 1.0e-8;
  check.gradMaxCoeff = 1.0e-7;
  check.gradRMS = 1.0e-8;
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
