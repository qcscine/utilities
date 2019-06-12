/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/GeometryOptimization/GeometryOptimizer.h"
#include "Utils/GeometryOptimization/IRCOptimizer.h"
#include "Utils/Optimizer/GradientBased/Bofill.h"
#include "Utils/Optimizer/GradientBased/LBFGS.h"
#include "Utils/Optimizer/GradientBased/SteepestDescent.h"
#include "Utils/Optimizer/HessianBased/EigenvectorFollowing.h"
#include "Utils/Optimizer/HessianBased/NewtonRaphson.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/CalculatorBasics/StatesHandler.h>
#include <Utils/Settings.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

// Define a mock states handler
class GeoOptMockStatesHandler : public StatesHandler {
 public:
  std::shared_ptr<State> getCurrentState(StateSize /*size*/) const override {
    return nullptr;
  };
};

// Define a mock calculator
class GeoOptMockCalculator : public Core::Calculator {
 public:
  GeoOptMockCalculator(GeoOptMockStatesHandler sh) : sh_(sh){};
  ~GeoOptMockCalculator() = default;
  void setStructure(const AtomCollection& structure) override {
    structure_ = structure;
  };
  void modifyPositions(PositionCollection newPositions) override {
    structure_.setPositions(newPositions);
  };
  const PositionCollection& getPositions() const override {
    return structure_.getPositions();
  };
  void setRequiredProperties(const PropertyList& requiredProperties) override{};
  PropertyList possibleProperties() const override {
    return Utils::Property::Energy | Utils::Property::Gradients;
  };
  const Results& calculate(std::string dummy = "") override {
    auto p1 = structure_.getPosition(0);
    auto p2 = structure_.getPosition(1);
    auto p3 = structure_.getPosition(2);
    auto v12 = p1 - p2;
    auto v32 = p3 - p2;
    auto r12 = v12.norm();
    auto r23 = v32.norm();

    // Energy
    double e = cos(r12) + cos(r23);
    // Gradient
    GradientCollection g(structure_.size(), 3);
    g(0, 0) = -sin(r12) * (v12[0] / r12);
    g(0, 1) = -sin(r12) * (v12[1] / r12);
    g(0, 2) = -sin(r12) * (v12[2] / r12);
    g(1, 0) = sin(r12) * (v12[0] / r12) + sin(r23) * (v32[0] / r23);
    g(1, 1) = sin(r12) * (v12[1] / r12) + sin(r23) * (v32[1] / r23);
    g(1, 2) = sin(r12) * (v12[2] / r12) + sin(r23) * (v32[2] / r23);
    g(2, 0) = -sin(r23) * (v32[0] / r23);
    g(2, 1) = -sin(r23) * (v32[1] / r23);
    g(2, 2) = -sin(r23) * (v32[2] / r23);

    // Hessian
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(9, 9);
    // diagonal
    h(0, 0) = (v12[0] * v12[0] * sin(r12) / (r12 * r12 * r12)) - (v12[0] * v12[0] * cos(r12) / (r12 * r12)) - (sin(r12) / r12);
    h(1, 1) = (v12[1] * v12[1] * sin(r12) / (r12 * r12 * r12)) - (v12[1] * v12[1] * cos(r12) / (r12 * r12)) - (sin(r12) / r12);
    h(2, 2) = (v12[2] * v12[2] * sin(r12) / (r12 * r12 * r12)) - (v12[2] * v12[2] * cos(r12) / (r12 * r12)) - (sin(r12) / r12);
    h(6, 6) = (v32[0] * v32[0] * sin(r23) / (r23 * r23 * r23)) - (v32[0] * v32[0] * cos(r23) / (r23 * r23)) - (sin(r23) / r23);
    h(7, 7) = (v32[1] * v32[1] * sin(r23) / (r23 * r23 * r23)) - (v32[1] * v32[1] * cos(r23) / (r23 * r23)) - (sin(r23) / r23);
    h(8, 8) = (v32[2] * v32[2] * sin(r23) / (r23 * r23 * r23)) - (v32[2] * v32[2] * cos(r23) / (r23 * r23)) - (sin(r23) / r23);
    h(3, 3) = h(0, 0) + h(6, 6);
    h(4, 4) = h(1, 1) + h(7, 7);
    h(5, 5) = h(2, 2) + h(8, 8);
    // rest of h-1-1
    h(1, 0) = ((sin(r12) / r12) - cos(r12)) * v12[1] * v12[0] / (r12 * r12);
    h(2, 0) = ((sin(r12) / r12) - cos(r12)) * v12[2] * v12[0] / (r12 * r12);
    h(2, 1) = ((sin(r12) / r12) - cos(r12)) * v12[2] * v12[1] / (r12 * r12);
    h(0, 1) = h(1, 0);
    h(0, 2) = h(2, 0);
    h(1, 2) = h(2, 1);
    // rest of h-3-3
    h(7, 6) = ((sin(r23) / r23) - cos(r23)) * v32[1] * v32[0] / (r23 * r23);
    h(8, 6) = ((sin(r23) / r23) - cos(r23)) * v32[2] * v32[0] / (r23 * r23);
    h(8, 7) = ((sin(r23) / r23) - cos(r23)) * v32[2] * v32[1] / (r23 * r23);
    h(6, 7) = h(7, 6);
    h(6, 8) = h(8, 6);
    h(7, 8) = h(8, 7);
    // rest of h-2-2
    h(4, 3) = h(1, 0) + h(7, 6);
    h(5, 3) = h(2, 0) + h(8, 6);
    h(5, 4) = h(2, 1) + h(8, 7);
    h(3, 4) = h(4, 3);
    h(3, 5) = h(5, 3);
    h(4, 5) = h(5, 4);
    // h-1-2 and h-2-1
    h(0, 3) = -h(0, 0);
    h(1, 4) = -h(1, 1);
    h(2, 5) = -h(2, 2);
    h(3, 0) = -h(0, 0);
    h(4, 1) = -h(1, 1);
    h(5, 2) = -h(2, 2);
    h(1, 3) = -h(1, 0);
    h(2, 3) = -h(2, 0);
    h(2, 4) = -h(2, 1);
    h(0, 4) = -h(0, 1);
    h(0, 5) = -h(0, 2);
    h(1, 5) = -h(1, 2);
    h(3, 1) = -h(1, 0);
    h(3, 2) = -h(2, 0);
    h(4, 2) = -h(2, 1);
    h(4, 0) = -h(0, 1);
    h(5, 0) = -h(0, 2);
    h(5, 1) = -h(1, 2);
    // h-3-2 and h-2-3
    h(6, 3) = -h(6, 6);
    h(7, 4) = -h(7, 7);
    h(8, 5) = -h(8, 8);
    h(3, 6) = -h(6, 6);
    h(4, 7) = -h(7, 7);
    h(5, 8) = -h(8, 8);
    h(7, 3) = -h(7, 6);
    h(8, 3) = -h(8, 6);
    h(8, 4) = -h(8, 7);
    h(6, 4) = -h(6, 7);
    h(6, 5) = -h(6, 8);
    h(7, 5) = -h(7, 8);
    h(3, 7) = -h(7, 6);
    h(3, 8) = -h(8, 6);
    h(4, 8) = -h(8, 7);
    h(4, 6) = -h(6, 7);
    h(5, 6) = -h(6, 8);
    h(5, 7) = -h(7, 8);

    r_ = Results();
    r_.setEnergy(e);
    r_.setGradients(g);
    r_.setHessian(h);
    return r_;
  };
  std::string name() const override {
    return std::string("GeoOptMockCalculator");
  };
  StatesHandler& statesHandler() override {
    return sh_;
  }
  const Settings& settings() const override {
    return settings_;
  }
  Settings& settings() override {
    return settings_;
  }
  const StatesHandler& statesHandler() const override {
    return sh_;
  }
  virtual Utils::Results& results() override {
    return r_;
  }
  virtual const Utils::Results& results() const override {
    return r_;
  }
  virtual std::unique_ptr<Utils::AtomCollection> getStructure() const override {
    return nullptr;
  }

 private:
  AtomCollection structure_;
  GeoOptMockStatesHandler sh_;
  Results r_;
  Settings settings_ = Settings("dummy");
  Core::Calculator* cloneImpl() const override {
    return nullptr;
  }
};

/**
 * @class Scine::Utils::Tests::GeometryOptimizerTests
 * @brief Comprises tests for the class Scine::Utils::GeometryOptimizer.
 * @test
 */
TEST(GeometryOptimizerTests, SteepestDescent) {
  GeoOptMockStatesHandler mockStatesHandler;
  GeoOptMockCalculator mockCalculator(mockStatesHandler);
  GeometryOptimizer<SteepestDescent> geo(mockCalculator);
  geo.check.maxIter = 300;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-9;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-9;
  geo.check.deltaValue = 1.0e-12;
  geo.check.requirement = 4;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -0.25 * M_PI;
  positions(2, 0) = 0.50 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = geo.optimize(atoms);

  // Check results
  EXPECT_TRUE(nIter < geo.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, M_PI, 1.0e-8);
  EXPECT_NEAR(r23, M_PI, 1.0e-8);
}

TEST(GeometryOptimizerTests, LBFGS) {
  GeoOptMockStatesHandler mockStatesHandler;
  GeoOptMockCalculator mockCalculator(mockStatesHandler);
  GeometryOptimizer<LBFGS> geo(mockCalculator);
  geo.check.maxIter = 50;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-9;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-9;
  geo.check.deltaValue = 1.0e-12;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -0.59 * M_PI;
  positions(2, 0) = 0.59 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = geo.optimize(atoms);

  // Check results
  EXPECT_TRUE(nIter < geo.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, M_PI, 1.0e-8);
  EXPECT_NEAR(r23, M_PI, 1.0e-8);
}

TEST(GeometryOptimizerTests, NewtonRaphson) {
  GeoOptMockStatesHandler mockStatesHandler;
  GeoOptMockCalculator mockCalculator(mockStatesHandler);
  GeometryOptimizer<NewtonRaphson> geo(mockCalculator);
  geo.check.maxIter = 100;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-8;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-8;
  geo.check.deltaValue = 1.0e-12;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -0.76 * M_PI;
  positions(2, 1) = 0.76 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = geo.optimize(atoms);

  // Check results
  EXPECT_TRUE(nIter < geo.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, M_PI, 1.0e-8);
  EXPECT_NEAR(r23, M_PI, 1.0e-8);
}

// TEST(GeometryOptimizerTests, Bofill) {
//   GeoOptMockStatesHandler mockStatesHandler;
//   GeoOptMockCalculator mockCalculator(mockStatesHandler);
//   GeometryOptimizer<Bofill> geo(mockCalculator);
//   geo.check.maxIter = 20000;
//   auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
//   Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
//   positions(0, 0) = -1.05 * M_PI;
//   positions(2, 1) = 0.85 * M_PI;
//   AtomCollection atoms(elements, positions);
//   auto nIter = geo.optimize(atoms);
//   // Check results
//   EXPECT_TRUE(nIter < geo.check.maxIter);
//   auto p1 = atoms.getPosition(0);
//   auto p2 = atoms.getPosition(1);
//   auto p3 = atoms.getPosition(2);
//   auto r12 = (p1 - p2).norm();
//   auto r23 = (p3 - p2).norm();
//   EXPECT_NEAR(r12, M_PI, 1.0e-8);
//   EXPECT_NEAR(r23, M_PI, 1.0e-8);
// }

// TEST(GeometryOptimizerTests, EigenvectorFollowing) {
//   GeoOptMockStatesHandler mockStatesHandler;
//   GeoOptMockCalculator mockCalculator(mockStatesHandler);
//   GeometryOptimizer<EigenvectorFollowing> geo(mockCalculator);
//   geo.check.maxIter = 20000;
//   auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
//   Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
//   positions(0, 0) = -0.75 * M_PI;
//   positions(2, 0) = 0.85 * M_PI;
//   AtomCollection atoms(elements, positions);
//   auto nIter = geo.optimize(atoms);
//   // Check results
//   EXPECT_TRUE(nIter < geo.check.maxIter);
//   auto p1 = atoms.getPosition(0);
//   auto p2 = atoms.getPosition(1);
//   auto p3 = atoms.getPosition(2);
//   auto r12 = (p1 - p2).norm();
//   auto r23 = (p3 - p2).norm();
//   EXPECT_NEAR(r12, M_PI, 1.0e-8);
//   EXPECT_NEAR(r23, M_PI, 1.0e-8);
// }

TEST(GeometryOptimizerTests, TestObserver) {
  GeoOptMockStatesHandler mockStatesHandler;
  GeoOptMockCalculator mockCalculator(mockStatesHandler);
  GeometryOptimizer<SteepestDescent> geo(mockCalculator);
  geo.check.maxIter = 200;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-8;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-8;
  geo.check.deltaValue = 1.0e-9;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -0.05 * M_PI;
  positions(2, 1) = 0.95 * M_PI;
  AtomCollection atoms(elements, positions);
  // Def observer
  int counter = 0;
  auto func = [&](const int&, const double&, const Eigen::VectorXd&) { counter++; };
  // Run
  geo.addObserver(func);
  auto nIter = geo.optimize(atoms);
  // Check
  EXPECT_EQ(nIter, counter);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
