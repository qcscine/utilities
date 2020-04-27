/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Optimizer/GradientBased/Dimer.h"
#include "Utils/Optimizer/GradientBased/Lbfgs.h"
#include "Utils/Optimizer/GradientBased/SteepestDescent.h"
#include "Utils/Optimizer/HessianBased/Bofill.h"
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

auto const optionalHessianTestFunction = [](const Eigen::VectorXd& parameters, double& value,
                                            Eigen::VectorXd& gradients, Eigen::MatrixXd& hessian, bool hessianUpdate) {
  value = -cos(parameters[0]) - 0.5 * cos(parameters[1]);
  gradients[0] = sin(parameters[0]);
  gradients[1] = 0.5 * sin(parameters[1]);
  if (hessianUpdate) {
    hessian(0, 0) = cos(parameters[0]);
    hessian(0, 1) = 0.0;
    hessian(1, 0) = 0.0;
    hessian(1, 1) = 0.5 * cos(parameters[1]);
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
  value = sin(parameters[0]) + sin(parameters[1]) + sin(parameters[2]) + sin(parameters[3]) + sin(parameters[4]);
  gradients[0] = cos(parameters[0]);
  gradients[1] = cos(parameters[1]);
  gradients[2] = cos(parameters[2]);
  gradients[3] = cos(parameters[3]);
  gradients[4] = cos(parameters[4]);
};

auto const sinSumHessianTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                          Eigen::MatrixXd& hessian, bool hessianUpdate) {
  value = sin(parameters[0]) + sin(parameters[1]) + sin(parameters[2]) + sin(parameters[3]) + sin(parameters[4]);
  gradients[0] = cos(parameters[0]);
  gradients[1] = cos(parameters[1]);
  gradients[2] = cos(parameters[2]);
  gradients[3] = cos(parameters[3]);
  gradients[4] = cos(parameters[4]);
  if (hessianUpdate) {
    hessian = Eigen::MatrixXd::Zero(5, 5);
    hessian(0, 0) = -sin(parameters[0]);
    hessian(1, 1) = -sin(parameters[1]);
    hessian(2, 2) = -sin(parameters[2]);
    hessian(3, 3) = -sin(parameters[3]);
    hessian(4, 4) = -sin(parameters[4]);
  }
};

double mullerBrownSingleTerm(double x, double y, double A, double a, double b, double c, double x0, double y0) {
  return A * exp(a * pow((x - x0), 2) + b * (x - x0) * (y - y0) + c * pow((y - y0), 2));
};

auto const mullerBrowmGradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
  value = 0.0;
  std::vector<double> Aconstants = {-200, -100, -170, 15};
  std::vector<double> aconstants = {-1, -1, -6.5, -0.7};
  std::vector<double> bconstants = {0, 0, 11, 0.6};
  std::vector<double> cconstants = {-10, -10, -6.5, 0.7};
  std::vector<double> x0constants = {1, 0, -0.5, -1};
  std::vector<double> y0constants = {0, 0.5, 1.5, 1};
  for (int i = 0; i < 4; ++i) {
    value += mullerBrownSingleTerm(parameters[0], parameters[1], Aconstants[i], aconstants[i], bconstants[i],
                                   cconstants[i], x0constants[i], y0constants[i]);
  }
};

auto const osdGradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
  value = pow((pow(parameters[0], 2) - 1), 2) + pow(parameters[1], 2);
  gradients[0] = 4 * parameters[0] * (pow(parameters[0], 2) - 1);
  gradients[1] = 2 * parameters[1];
};

auto const osdHessianTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                       Eigen::MatrixXd& hessian, bool hessianUpdate) {
  value = pow((pow(parameters[0], 2) - 1), 2) + pow(parameters[1], 2);
  gradients[0] = 4 * parameters[0] * (pow(parameters[0], 2) - 1);
  gradients[1] = 2 * parameters[1];
  if (hessianUpdate) {
    hessian(0, 0) = 12 * pow(parameters[0], 2) - 4;
    hessian(1, 1) = 2;
    hessian(0, 1) = 0;
    hessian(1, 0) = 0;
  }
};

auto const minyaevQuappGradientTestFunction = [](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
  value = cos(2 * parameters[0]) + cos(2 * parameters[1]) + 0.57 * cos(2 * parameters[0] - 2 * parameters[1]);
  gradients[0] = -1.14 * sin(2 * parameters[0] - 2 * parameters[1]) - 2 * sin(2 * parameters[0]);
  gradients[1] = 1.14 * sin(2 * parameters[0] - 2 * parameters[1]) - 2 * sin(2 * parameters[1]);
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
TEST(OptimizerTests, LbfgsTest) {
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check);
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check);
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check);
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check);
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check);
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, EigenvectorFollowingTest_4) {
  EigenvectorFollowing optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.7;
  positions[1] = -0.5;
  GradientBasedCheck check;
  check.maxIter = 500;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, shang19HessianTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, EigenvectorFollowingTest_5) {
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
  auto nCycles = optimizer.optimize(positions, sinSumHessianTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[1], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[2], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[3], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[4], M_PI / 2, 1.0e-8);
}
TEST(OptimizerTests, EigenvectorFollowingTest_6) {
  EigenvectorFollowing optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.6;
  positions[1] = 1;
  GradientBasedCheck check;
  check.maxIter = 500;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, osdHessianTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}
TEST(OptimizerTests, DimerTest_1) {
  Dimer optimizer;
  optimizer.useGdiis = true;
  Eigen::VectorXd positions(2);
  positions[0] = 0.49 * M_PI;
  positions[1] = 0.51 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-5);
  EXPECT_NEAR(positions[1], M_PI, 1.0e-5);
}
TEST(OptimizerTests, DimerTest_2) {
  Dimer optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.51 * M_PI;
  positions[1] = 0.49 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-6);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-6);
}
TEST(OptimizerTests, DimerTest_3) {
  Dimer optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.51 * M_PI;
  positions[1] = 0.51 * M_PI;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, gradientTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-5);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-5);
}
TEST(OptimizerTests, DimerTest_4) {
  Dimer optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.7;
  positions[1] = -0.5;
  GradientBasedCheck check;
  check.maxIter = 400;
  check.deltaValue = 1e-12;
  optimizer.defaultTranslationStep = 0.01;
  auto nCycles = optimizer.optimize(positions, shang19GradientTestFunction, check);
  EXPECT_TRUE(nCycles < 400);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-5);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-5);
}
TEST(OptimizerTests, DimerTest_5) {
  Dimer optimizer;
  Eigen::VectorXd positions(5);
  positions[0] = -0.1;
  positions[1] = -0.2;
  positions[2] = -0.3;
  positions[3] = -0.4;
  positions[4] = 0.5;
  GradientBasedCheck check;
  check.maxIter = 50;
  check.deltaValue = 1e-12;
  auto nCycles = optimizer.optimize(positions, sinSumGradientTestFunction, check);
  EXPECT_TRUE(nCycles < 50);
  EXPECT_NEAR(positions[0], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[1], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[2], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[3], -M_PI / 2, 1.0e-8);
  EXPECT_NEAR(positions[4], M_PI / 2, 1.0e-8);
}
TEST(OptimizerTests, DimerTest_6) {
  Dimer optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 0.0;
  positions[1] = 0.5;
  GradientBasedCheck check;
  check.maxIter = 100;
  check.deltaValue = 1e-12;
  optimizer.defaultTranslationStep = 0.1;
  auto nCycles = optimizer.optimize(positions, osdGradientTestFunction, check);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[0], 0.0, 1.0e-5);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-5);
}
TEST(OptimizerTests, DimerTest_7) {
  Dimer optimizer;
  Eigen::VectorXd positions(2);
  positions[0] = 1.4;
  positions[1] = 1.4;
  GradientBasedCheck check;
  check.maxIter = 400;
  check.deltaValue = 1e-12;
  optimizer.defaultTranslationStep = 0.1;
  auto nCycles = optimizer.optimize(positions, minyaevQuappGradientTestFunction, check);
  EXPECT_TRUE(nCycles < 400);
  EXPECT_NEAR(positions[0], M_PI / 2, 1.0e-5);
  EXPECT_NEAR(positions[1], M_PI / 2, 1.0e-5);
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
}
TEST(OptimizerTests, BofillTest_2) {
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check);
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
  auto nCycles = optimizer.optimize(positions, optionalHessianTestFunction, check);
  EXPECT_TRUE(nCycles < 100);
  EXPECT_NEAR(positions[0], M_PI, 1.0e-8);
  EXPECT_NEAR(positions[1], 0.0, 1.0e-8);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
