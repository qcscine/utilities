/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Optimizer/GradientBased/Bfgs.h"
#include "Utils/Optimizer/GradientBased/Dimer.h"
#include "Utils/Optimizer/GradientBased/Lbfgs.h"
#include "Utils/Optimizer/GradientBased/SteepestDescent.h"
#include "Utils/Optimizer/HessianBased/Bofill.h"
#include "Utils/Optimizer/HessianBased/EigenvectorFollowing.h"
#include "Utils/Optimizer/HessianBased/NewtonRaphson.h"
#include <Core/Interfaces/Calculator.h>
#include <Core/Log.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Settings.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

auto const gradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
  value = -std::cos(parameters[0]) - 0.5 * std::cos(parameters[1]);
  gradients[0] = std::sin(parameters[0]);
  gradients[1] = 0.5 * std::sin(parameters[1]);
};

auto const optionalHessianTestFunction = [](const Eigen::VectorXd& parameters, double& value,
                                            Eigen::VectorXd& gradients, Eigen::MatrixXd& hessian, bool hessianUpdate) {
  value = -std::cos(parameters[0]) - 0.5 * std::cos(parameters[1]);
  gradients[0] = std::sin(parameters[0]);
  gradients[1] = 0.5 * std::sin(parameters[1]);
  if (hessianUpdate) {
    hessian(0, 0) = std::cos(parameters[0]);
    hessian(0, 1) = 0.0;
    hessian(1, 0) = 0.0;
    hessian(1, 1) = 0.5 * std::cos(parameters[1]);
  }
};

auto const shang19GradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
  value = pow(parameters[0], 4) + 4 * pow(parameters[0], 2) * pow(parameters[1], 2) - 2 * (parameters[0], 2) +
          2 * pow(parameters[1], 2);
  gradients[0] = 4 * parameters[0] * (pow(parameters[0], 2) + 2 * pow(parameters[1], 2) - 1);
  gradients[1] = 4 * (2 * pow(parameters[0], 2) + 1) * parameters[1];
};

auto const shang19HessianTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                           Eigen::MatrixXd& hessian, bool hessianUpdate) {
  value = pow(parameters[0], 4) + 4 * pow(parameters[0], 2) * pow(parameters[1], 2) - 2 * (parameters[0], 2) +
          2 * pow(parameters[1], 2);
  gradients[0] = 4 * parameters[0] * (pow(parameters[0], 2) + 2 * pow(parameters[1], 2) - 1);
  gradients[1] = 4 * (2 * pow(parameters[0], 2) + 1) * parameters[1];
  if (hessianUpdate) {
    hessian(0, 0) = 12 * pow(parameters[0], 2) + 8 * pow(parameters[1], 2) - 4;
    hessian(0, 1) = 16 * parameters[0] * parameters[1];
    hessian(1, 0) = 16 * parameters[0] * parameters[1];
    hessian(1, 1) = 8 * pow(parameters[0], 2) + 4;
  }
};

auto const sinSumGradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
  value = std::sin(parameters[0]) + std::sin(parameters[1]) + std::sin(parameters[2]) + std::sin(parameters[3]) +
          std::sin(parameters[4]);
  gradients[0] = std::cos(parameters[0]);
  gradients[1] = std::cos(parameters[1]);
  gradients[2] = std::cos(parameters[2]);
  gradients[3] = std::cos(parameters[3]);
  gradients[4] = std::cos(parameters[4]);
};

auto const sinSumHessianTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                          Eigen::MatrixXd& hessian, bool hessianUpdate) {
  value = std::sin(parameters[0]) + std::sin(parameters[1]) + std::sin(parameters[2]) + std::sin(parameters[3]) +
          std::sin(parameters[4]);
  gradients[0] = std::cos(parameters[0]);
  gradients[1] = std::cos(parameters[1]);
  gradients[2] = std::cos(parameters[2]);
  gradients[3] = std::cos(parameters[3]);
  gradients[4] = std::cos(parameters[4]);
  if (hessianUpdate) {
    hessian = Eigen::MatrixXd::Zero(5, 5);
    hessian(0, 0) = -std::sin(parameters[0]);
    hessian(1, 1) = -std::sin(parameters[1]);
    hessian(2, 2) = -std::sin(parameters[2]);
    hessian(3, 3) = -std::sin(parameters[3]);
    hessian(4, 4) = -std::sin(parameters[4]);
  }
};

auto const osdGradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
  value = pow((pow(parameters[0], 2) - 1), 2) + pow(parameters[1], 2);
  gradients[0] = 4 * parameters[0] * (pow(parameters[0], 2) - 1);
  gradients[1] = 2 * parameters[1];
};

auto const osdHessianTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                       Eigen::MatrixXd& hessian, bool hessianUpdate) {
  value = std::pow((std::pow(parameters[0], 2) - 1), 2) + std::pow(parameters[1], 2);
  gradients[0] = 4 * parameters[0] * (std::pow(parameters[0], 2) - 1);
  gradients[1] = 2 * parameters[1];
  if (hessianUpdate) {
    hessian(0, 0) = 12 * std::pow(parameters[0], 2) - 4;
    hessian(1, 1) = 2;
    hessian(0, 1) = 0;
    hessian(1, 0) = 0;
  }
};

auto const minyaevQuappGradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
  value = std::cos(2 * parameters[0]) + std::cos(2 * parameters[1]) + 0.57 * std::cos(2 * parameters[0] - 2 * parameters[1]);
  gradients[0] = -1.14 * std::sin(2 * parameters[0] - 2 * parameters[1]) - 2 * std::sin(2 * parameters[0]);
  gradients[1] = 1.14 * std::sin(2 * parameters[0] - 2 * parameters[1]) - 2 * std::sin(2 * parameters[1]);
};

/**
 * @class Scine::Utils::Tests::OptimizerTests
 * @brief Comprises tests for the class Scine::Utils::Optimizer derived classes.
 * @test
 */
TEST(OptimizerTests, SDTest) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
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
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check, log);

  EXPECT_TRUE(nCycles < 1000);
  EXPECT_NEAR(positions[0], 0.0, 2.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 2.0e-8);
}
TEST(OptimizerTests, DynamicSDTest) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  SteepestDescent optimizer;
  optimizer.useTrustRadius = true;
  optimizer.trustRadius = 0.05;
  optimizer.dynamicMultiplier = 1.1;
  Eigen::VectorXd positions(2);
  positions[0] = 0.25 * M_PI;
  positions[1] = 0.75 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 100;
  check.deltaValue = 1e-12;
  check.stepMaxCoeff = 1.0e-7;
  check.stepRMS = 1.0e-9;
  check.gradMaxCoeff = 1.0e-7;
  check.gradRMS = 1.0e-9;
  check.requirement = 4;
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check, log);

  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[0], 0.0, 2.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 2.0e-8);
}
TEST(OptimizerTests, LbfgsTest) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  Lbfgs optimizer;
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
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check, log);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, NewtonRaphsonTest_Minimization) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check, log);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, NewtonRaphsonTest_PartialMaximization) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check, log);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], M_PI, 1.0e-8);
}
TEST(OptimizerTests, NewtonRaphsonTest_Maximization) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  NewtonRaphson optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.51 * M_PI;
  positions[1] = 0.51 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check, log);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], M_PI, 1.0e-8);
}
TEST(OptimizerTests, EigenvectorFollowingTest_1) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  EigenvectorFollowing optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.49 * M_PI;
  positions[1] = 0.51 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check, log);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], M_PI, 1.0e-8);
}
TEST(OptimizerTests, EigenvectorFollowingTest_2) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  EigenvectorFollowing optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.51 * M_PI;
  positions[1] = 0.49 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check, log);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, EigenvectorFollowingTest_3) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  EigenvectorFollowing optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.51 * M_PI;
  positions[1] = 0.51 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check, log);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, EigenvectorFollowingTest_4) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  EigenvectorFollowing optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.7;
  positions[1] = -0.5;
  GradientBasedCheck check;
  check.maxIter = 500;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, shang19HessianTestFunction, check, log);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, EigenvectorFollowingTest_5) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  EigenvectorFollowing optimizer;
  Eigen::VectorXd positions(5);
  positions[0] = -0.1;
  positions[1] = -0.2;
  positions[2] = -0.3;
  positions[3] = -0.4;
  positions[4] = 0.5;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, sinSumHessianTestFunction, check, log);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[1], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[2], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[3], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[4], M_PI / 2, 1.0e-8);
}
TEST(OptimizerTests, EigenvectorFollowingTest_6) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  EigenvectorFollowing optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.6;
  positions[1] = 1;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, osdHessianTestFunction, check, log);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, DimerTest_1) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  Dimer optimizer;
  Eigen::VectorXd positions(2);
  optimizer.translationMethod = "AMSGRAD";
  positions[0] = 0.49 * M_PI;
  positions[1] = 0.51 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 550;
  check.deltaValue = 1e-14;
  check.gradMaxCoeff = 1e-10;
  check.requirement = 4;
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check, log);
  EXPECT_TRUE(nCycles < 550);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], M_PI, 1.0e-8);
}
TEST(OptimizerTests, DimerTest_2) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  Dimer optimizer;
  optimizer.translationMethod = "Linesearch";
  Eigen::VectorXd positions(2);
  positions[0] = 0.51 * M_PI;
  positions[1] = 0.49 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 100;
  check.deltaValue = 1e-14;
  check.gradMaxCoeff = 1e-10;
  check.requirement = 4;
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check, log);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, DimerTest_3) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  Dimer optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.51 * M_PI;
  positions[1] = 0.51 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check, log);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, DimerTest_4) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  Dimer optimizer;
  optimizer.maxValueMemory = 10;
  Eigen::VectorXd positions(2);
  positions[0] = 0.7;
  positions[1] = -0.5;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  optimizer.defaultTranslationStep = 1.0;
  auto nCycles = optimizer.optimize(positions, shang19GradientTestFunction, check, log);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, DimerTest_5) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  Dimer optimizer;
  optimizer.translationMethod = "AMSGRAD";
  Eigen::VectorXd positions(5);
  positions[0] = -0.1;
  positions[1] = -0.2;
  positions[2] = -0.3;
  positions[3] = -0.4;
  positions[4] = 0.5;
  GradientBasedCheck check;
  check.maxIter = 500;
  check.deltaValue = 1e-14;
  check.gradRMS = 1.0e-10;
  check.requirement = 4;
  auto nCycles = optimizer.optimize(positions, sinSumGradientTestFunction, check, log);
  EXPECT_TRUE(nCycles < 500);
  EXPECT_NEAR(positions[0], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[1], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[2], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[3], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[4], M_PI / 2, 1.0e-8);
}
TEST(OptimizerTests, DimerTest_6) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  Dimer optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.0;
  positions[1] = 0.5;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-14;
  check.gradMaxCoeff = 1e-10;
  check.requirement = 4;
  optimizer.defaultTranslationStep = 0.1;
  auto nCycles = optimizer.optimize(positions, osdGradientTestFunction, check, log);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, DimerTest_7) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  Dimer optimizer;
  Eigen::VectorXd positions(2);
  optimizer.skipFirstRotation = false;
  optimizer.decreaseRotationGradientThreshold = true;
  optimizer.gradientInterpolation = false;
  optimizer.rotationCG = false;
  optimizer.onlyOneRotation = false;
  optimizer.translationMethod = "BFGS";
  optimizer.multiScale = false;
  optimizer.radius = 0.02;
  optimizer.phiTolerance = 1e-2;
  optimizer.rotationGradientThresholdFirstCycle = 2e-7;
  optimizer.rotationGradientThresholdOtherCycles = 2e-4;
  optimizer.loweredRotationGradientThreshold = 2e-3;
  optimizer.gradRMSDthreshold = 1e-2;
  optimizer.trustRadius = 0.1;
  optimizer.defaultTranslationStep = 0.1;
  optimizer.alpha1 = 0.02;
  optimizer.beta1 = 0.2;
  optimizer.beta2 = 0.02;
  optimizer.maxRotationsFirstCycle = 200;
  optimizer.maxRotationsOtherCycles = 200;
  optimizer.intervalOfRotations = 4;
  optimizer.cycleOfRotationGradientDecrease = 4;
  optimizer.lbfgsMemory = 7;
  optimizer.bfgsStart = 10;
  optimizer.minimizationCycle = 10;
  positions[0] = 1.4;
  positions[1] = 1.4;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  check.gradMaxCoeff = 1e-10;
  check.requirement = 4;
  auto nCycles = optimizer.optimize(positions, minyaevQuappGradientTestFunction, check, log);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[1], M_PI / 2, 1.0e-8);
}
TEST(OptimizerTests, BofillTest_1) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check, log);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
}
TEST(OptimizerTests, BofillTest_2) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  Bofill optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.51 * M_PI;
  positions[1] = 0.49 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 100;
  check.deltaValue = 1e-12;
  check.stepMaxCoeff = 1.0e-7;
  check.stepRMS = 1.0e-8;
  check.gradMaxCoeff = 1.0e-7;
  check.gradRMS = 1.0e-8;
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check, log);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, BofillTest_3) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check, log);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
}
TEST(OptimizerTests, BofillTest_4) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check, log);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, OscillationTest) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  Bofill optimizer;
  Eigen::VectorXd positions1(2);
  Eigen::VectorXd positions2(2);
  positions1[0] = 0.50;
  positions1[1] = 0.50;
  positions2[0] = -1 * positions1[0];
  positions2[1] = positions1[1];
  optimizer.maxValueMemory = 3;
  EXPECT_FALSE(optimizer.isOscillating(4.0));
  EXPECT_FALSE(optimizer.isOscillating(2.0));
  EXPECT_TRUE(optimizer.isOscillating(4.0));
  optimizer.oscillationCorrection((positions2 - positions1), positions2);
  EXPECT_NEAR(positions2[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions2[1], 0.5, 1.0e-8);
}

TEST(OptimizerTests, PrepareRestartTest_SD) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  SteepestDescent optimizer;
  // Add an observer storing the cycle numbers
  std::vector<int> cycleNumbers;
  auto func = [&](const int& cycle, const double& /*energy*/, const Eigen::VectorXd& /*params*/) {
    cycleNumbers.push_back(cycle);
  };
  optimizer.addObserver(func);

  Eigen::VectorXd positions(2);
  positions[0] = 0.25 * M_PI;
  positions[1] = 0.75 * M_PI;
  GradientBasedCheck check;

  // Set start cycle to 4 and finish after sixth cycle
  optimizer.prepareRestart(4);
  check.maxIter = 6;

  optimizer.optimize(positions, gradientTestFunction, check, log);
  ASSERT_THAT(cycleNumbers, ElementsAre(4, 5, 6));
}

TEST(OptimizerTests, PrepareRestartTest_Bfgs) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  Bfgs optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.25 * M_PI;
  positions[1] = 0.75 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 2;
  optimizer.optimize(positions, gradientTestFunction, check, log);
  // Prepare restart
  optimizer.prepareRestart(3);
  // Inverse Hessian should have been cleared
  EXPECT_TRUE(optimizer.invH.size() == 0);
  // Add an observer storing the cycle numbers
  std::vector<int> cycleNumbers;
  auto func = [&](const int& cycle, const double& /*energy*/, const Eigen::VectorXd& /*params*/) {
    cycleNumbers.push_back(cycle);
  };
  optimizer.addObserver(func);
  // Restart and check cycle counts
  check.maxIter = 4;
  optimizer.optimize(positions, gradientTestFunction, check, log);
  ASSERT_THAT(cycleNumbers, ElementsAre(3, 4));
}

TEST(OptimizerTests, PrepareRestartTest_Dimer) {
  auto log = Core::Log::silent();
  log.output.add("cout", Core::Log::coutSink());
  log.warning.add("cerr", Core::Log::cerrSink());
  log.error.add("cerr", Core::Log::cerrSink());
  Dimer optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.25 * M_PI;
  positions[1] = 0.75 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 2;
  optimizer.optimize(positions, gradientTestFunction, check, log);
  // Prepare restart
  optimizer.prepareRestart(1);
  // Inverse Hessian should have been cleared
  EXPECT_TRUE(optimizer.invH.size() == 0);
  // Dimer guess vector should be cleared
  EXPECT_FALSE(optimizer.guessVector);
  // Add an observer storing the cycle numbers
  std::vector<int> cycleNumbers;
  auto func = [&](const int& cycle, const double& /*energy*/, const Eigen::VectorXd& /*params*/) {
    cycleNumbers.push_back(cycle);
  };
  optimizer.addObserver(func);
  // Restart and check cycles
  optimizer.optimize(positions, gradientTestFunction, check, log);
  ASSERT_THAT(cycleNumbers, ElementsAre(1, 2));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
