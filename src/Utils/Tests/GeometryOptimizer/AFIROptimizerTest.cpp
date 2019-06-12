/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/GeometryOptimization/AFIROptimizer.h"
#include "Utils/Optimizer/GradientBased/LBFGS.h"
#include "Utils/Optimizer/GradientBased/SteepestDescent.h"
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
class AFIROptMockStatesHandler : public StatesHandler {
 public:
  std::shared_ptr<State> getCurrentState(StateSize /*size*/) const override {
    return nullptr;
  };
};

// Define a mock calculator
class AFIROptMockCalculator : public Core::Calculator {
 public:
  AFIROptMockCalculator(AFIROptMockStatesHandler sh) : sh_(sh){};
  ~AFIROptMockCalculator() = default;
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

    r_ = Results();
    r_.setEnergy(e);
    r_.setGradients(g);
    return r_;
  };
  std::string name() const override {
    return std::string("AFIROptMockCalculator");
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
    return std::make_unique<AtomCollection>(structure_);
  }

 private:
  AtomCollection structure_;
  AFIROptMockStatesHandler sh_;
  Results r_;
  Settings settings_ = Settings("dummy");
  Core::Calculator* cloneImpl() const override {
    return nullptr;
  }
};

/**
 * @class Scine::Utils::Tests::AFIROptimizerTests
 * @brief Comprises tests for the class Scine::Utils::AFIROptimizer.
 * @test
 */
TEST(AFIROptimizerTests, SteepestDescent) {
  AFIROptMockStatesHandler mockStatesHandler;
  AFIROptMockCalculator mockCalculator(mockStatesHandler);
  AFIROptimizer<SteepestDescent> afir(mockCalculator);
  afir.check.maxIter = 300;
  afir.check.stepMaxCoeff = 1.0e-7;
  afir.check.stepRMS = 1.0e-8;
  afir.check.gradMaxCoeff = 1.0e-7;
  afir.check.gradRMS = 1.0e-8;
  afir.check.deltaValue = 1.0e-12;
  afir.check.requirement = 4;
  afir.phaseIn = 0;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -0.25 * M_PI;
  positions(2, 0) = 0.50 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = afir.optimize(atoms);

  // Check results
  EXPECT_TRUE(nIter < afir.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, M_PI, 2.0e-8);
  EXPECT_NEAR(r23, M_PI, 2.0e-8);
}

TEST(AFIROptimizerTests, LBFGS_NoPotential) {
  AFIROptMockStatesHandler mockStatesHandler;
  AFIROptMockCalculator mockCalculator(mockStatesHandler);
  AFIROptimizer<LBFGS> afir(mockCalculator);
  afir.check.maxIter = 50;
  afir.check.stepMaxCoeff = 1.0e-7;
  afir.check.stepRMS = 1.0e-8;
  afir.check.gradMaxCoeff = 1.0e-7;
  afir.check.gradRMS = 1.0e-8;
  afir.check.deltaValue = 1.0e-9;
  afir.phaseIn = 0;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -0.59 * M_PI;
  positions(2, 0) = 0.59 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = afir.optimize(atoms);

  // Check results
  EXPECT_TRUE(nIter < afir.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, M_PI, 1.0e-8);
  EXPECT_NEAR(r23, M_PI, 1.0e-8);
}

TEST(AFIROptimizerTests, LBFGS_WeakOnly) {
  AFIROptMockStatesHandler mockStatesHandler;
  AFIROptMockCalculator mockCalculator(mockStatesHandler);
  AFIROptimizer<LBFGS> afir(mockCalculator);
  afir.check.maxIter = 50;
  afir.check.stepMaxCoeff = 1.0e-7;
  afir.check.stepRMS = 1.0e-8;
  afir.check.gradMaxCoeff = 1.0e-7;
  afir.check.gradRMS = 1.0e-8;
  afir.check.deltaValue = 1.0e-9;
  afir.phaseIn = 0;
  afir.weak = true;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -0.59 * M_PI;
  positions(2, 0) = 0.59 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = afir.optimize(atoms);

  // Check results
  EXPECT_TRUE(nIter < afir.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_TRUE(r12 < M_PI);
  EXPECT_TRUE(r23 < M_PI);
}

TEST(AFIROptimizerTests, LBFGS_ListAttractiveOnly) {
  AFIROptMockStatesHandler mockStatesHandler;
  AFIROptMockCalculator mockCalculator(mockStatesHandler);
  AFIROptimizer<LBFGS> afir(mockCalculator);
  afir.check.maxIter = 50;
  afir.check.stepMaxCoeff = 1.0e-7;
  afir.check.stepRMS = 1.0e-8;
  afir.check.gradMaxCoeff = 1.0e-7;
  afir.check.gradRMS = 1.0e-8;
  afir.check.deltaValue = 1.0e-9;
  afir.phaseIn = 0;
  afir.rhsList = {0};
  afir.lhsList = {2};
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -0.59 * M_PI;
  positions(2, 0) = 0.59 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = afir.optimize(atoms);

  // Check results
  EXPECT_TRUE(nIter < afir.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, 2.9814702666565984, 1.0e-8);
  EXPECT_NEAR(r23, 2.9814702666565984, 1.0e-8);
}

TEST(AFIROptimizerTests, LBFGS_ListRepulsiveOnly) {
  AFIROptMockStatesHandler mockStatesHandler;
  AFIROptMockCalculator mockCalculator(mockStatesHandler);
  AFIROptimizer<LBFGS> afir(mockCalculator);
  afir.check.maxIter = 50;
  afir.check.stepMaxCoeff = 1.0e-7;
  afir.check.stepRMS = 1.0e-8;
  afir.check.gradMaxCoeff = 1.0e-7;
  afir.check.gradRMS = 1.0e-8;
  afir.check.deltaValue = 1.0e-9;
  afir.attractive = false;
  afir.phaseIn = 0;
  afir.rhsList = {0};
  afir.lhsList = {2};
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -0.59 * M_PI;
  positions(2, 0) = 0.59 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = afir.optimize(atoms);

  // Check results
  EXPECT_TRUE(nIter < afir.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, 3.3017150393650367, 1.0e-8);
  EXPECT_NEAR(r23, 3.3017150393650367, 1.0e-8);
}

TEST(AFIROptimizerTests, LBFGS_ListAndWeak) {
  AFIROptMockStatesHandler mockStatesHandler;
  AFIROptMockCalculator mockCalculator(mockStatesHandler);
  AFIROptimizer<LBFGS> afir(mockCalculator);
  afir.check.maxIter = 50;
  afir.check.stepMaxCoeff = 1.0e-7;
  afir.check.stepRMS = 1.0e-8;
  afir.check.gradMaxCoeff = 1.0e-7;
  afir.check.gradRMS = 1.0e-8;
  afir.check.deltaValue = 1.0e-9;
  afir.phaseIn = 0;
  afir.weak = true;
  afir.transformCoordinates = false;
  afir.rhsList = {0};
  afir.lhsList = {2};
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -0.59 * M_PI;
  positions(2, 0) = 0.59 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = afir.optimize(atoms);

  // Check results
  EXPECT_TRUE(nIter < afir.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_TRUE(r12 < 2.9814702666565984);
  EXPECT_TRUE(r23 < 2.9814702666565984);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
