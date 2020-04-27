/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/GeometryOptimization/AfirOptimizer.h"
#include "Utils/Optimizer/GradientBased/Lbfgs.h"
#include "Utils/Optimizer/GradientBased/SteepestDescent.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/CalculatorBasics/StatesHandler.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

// Define a mock calculator
class AfirOptMockCalculator : public CloneInterface<AfirOptMockCalculator, Core::Calculator> {
 public:
  AfirOptMockCalculator() = default;
  ~AfirOptMockCalculator() final = default;
  void setStructure(const AtomCollection& structure) final {
    structure_ = structure;
  };
  void modifyPositions(PositionCollection newPositions) final {
    structure_.setPositions(newPositions);
  };
  const PositionCollection& getPositions() const final {
    return structure_.getPositions();
  };
  void setRequiredProperties(const PropertyList& requiredProperties) final{};
  PropertyList getRequiredProperties() const final {
    return Utils::PropertyList{};
  };
  PropertyList possibleProperties() const final {
    return Utils::Property::Energy | Utils::Property::Gradients;
  };
  const Results& calculate(std::string dummy = "") final {
    auto p1 = structure_.getPosition(0);
    auto p2 = structure_.getPosition(1);
    auto p3 = structure_.getPosition(2);
    auto v12 = p1 - p2;
    auto v32 = p3 - p2;
    auto r12 = v12.norm();
    auto r23 = v32.norm();

    // Energy
    // f(r12, r23) = (r12-4)*(r12-2)*(r12+2)*(r12+4) + (r23-4)*(r23-2)*(r23+2)*(r23+4);
    double e = pow(r12, 4) - 20 * pow(r12, 2) + pow(r23, 4) - 20 * pow(r23, 2) + 128;
    // Gradient
    GradientCollection g(structure_.size(), 3);
    g(0, 0) = (4.0 * pow(r12, 2.0) - 40.0) * v12[0];
    g(0, 1) = (4.0 * pow(r12, 2.0) - 40.0) * v12[1];
    g(0, 2) = (4.0 * pow(r12, 2.0) - 40.0) * v12[2];
    g(1, 0) = -(4.0 * pow(r12, 2.0) - 40.0) * v12[0] - (4.0 * pow(r23, 2.0) - 40.0) * v32[0];
    g(1, 1) = -(4.0 * pow(r12, 2.0) - 40.0) * v12[1] - (4.0 * pow(r23, 2.0) - 40.0) * v32[1];
    g(1, 2) = -(4.0 * pow(r12, 2.0) - 40.0) * v12[2] - (4.0 * pow(r23, 2.0) - 40.0) * v32[2];
    g(2, 0) = (4.0 * pow(r23, 2.0) - 40.0) * v32[0];
    g(2, 1) = (4.0 * pow(r23, 2.0) - 40.0) * v32[1];
    g(2, 2) = (4.0 * pow(r23, 2.0) - 40.0) * v32[2];

    r_ = Results();
    r_.set<Property::Energy>(e);
    r_.set<Property::Gradients>(g);
    return r_;
  };
  std::string name() const final {
    return std::string("AfirOptMockCalculator");
  };
  const Settings& settings() const final {
    return settings_;
  }
  Settings& settings() final {
    return settings_;
  }
  std::shared_ptr<Core::State> getState() const final {
    return nullptr;
  }
  void loadState(std::shared_ptr<Core::State> state) final {
  }
  Utils::Results& results() final {
    return r_;
  }
  const Utils::Results& results() const final {
    return r_;
  }
  std::unique_ptr<Utils::AtomCollection> getStructure() const final {
    return std::make_unique<AtomCollection>(structure_);
  }
  bool supportsMethodFamily(const std::string& methodFamily) const final {
    return true;
  }

 private:
  AtomCollection structure_;
  Results r_;
  Settings settings_ = Settings("dummy");
};

/**
 * @class Scine::Utils::Tests::AfirOptimizerTests
 * @brief Comprises tests for the class Scine::Utils::AfirOptimizer.
 * @test
 */
TEST(AfirOptimizerTests, SteepestDescent_NoPotential) {
  AfirOptMockCalculator mockCalculator;
  AfirOptimizer<SteepestDescent> afir(mockCalculator);
  afir.check.maxIter = 50;
  afir.optimizer.factor = 0.01;
  afir.check.stepMaxCoeff = 1.0e-7;
  afir.check.stepRMS = 1.0e-9;
  afir.check.gradMaxCoeff = 1.0e-7;
  afir.check.gradRMS = 1.0e-9;
  afir.check.deltaValue = 1.0e-12;
  afir.check.requirement = 4;
  afir.phaseIn = 0;
  afir.weak = false;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -M_PI;
  positions(2, 0) = M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = afir.optimize(atoms);

  // Check results
  EXPECT_TRUE(nIter < afir.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, sqrt(10), 2.0e-8);
  EXPECT_NEAR(r23, sqrt(10), 2.0e-8);
}

TEST(AfirOptimizerTests, Lbfgs_NoPotential) {
  AfirOptMockCalculator mockCalculator;
  AfirOptimizer<Lbfgs> afir(mockCalculator);
  afir.optimizer.maxBacktracking = 20;
  afir.transformCoordinates = false;
  afir.check.maxIter = 300;
  afir.check.stepMaxCoeff = 1.0e-7;
  afir.check.stepRMS = 1.0e-8;
  afir.check.gradMaxCoeff = 1.0e-7;
  afir.check.gradRMS = 1.0e-8;
  afir.check.deltaValue = 1.0e-12;
  afir.phaseIn = 0;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = 0.75 * M_PI;
  positions(2, 0) = 0.64 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = afir.optimize(atoms);

  // Check results
  EXPECT_TRUE(nIter < afir.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, sqrt(10), 1.0e-8);
  EXPECT_NEAR(r23, sqrt(10), 1.0e-8);
}

TEST(AfirOptimizerTests, Lbfgs_WeakOnly) {
  AfirOptMockCalculator mockCalculator;
  AfirOptimizer<Lbfgs> afir(mockCalculator);
  afir.check.maxIter = 50;
  afir.optimizer.linesearch = false;
  afir.transformCoordinates = false;
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
  EXPECT_TRUE(r12 < sqrt(10));
  EXPECT_TRUE(r23 < sqrt(10));
}

TEST(AfirOptimizerTests, Lbfgs_ListAttractiveOnly) {
  AfirOptMockCalculator mockCalculator;
  AfirOptimizer<Lbfgs> afir(mockCalculator);
  afir.check.maxIter = 50;
  afir.optimizer.linesearch = false;
  afir.transformCoordinates = false;
  afir.check.stepMaxCoeff = 1.0e-7;
  afir.check.stepRMS = 1.0e-8;
  afir.check.gradMaxCoeff = 1.0e-7;
  afir.check.gradRMS = 1.0e-8;
  afir.check.deltaValue = 1.0e-9;
  afir.phaseIn = 0;
  afir.rhsList = {0};
  afir.lhsList = {1};
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
  EXPECT_TRUE(r12 < sqrt(10));
  EXPECT_NEAR(r23, sqrt(10), 1.0e-8);
}

TEST(AfirOptimizerTests, Lbfgs_ListRepulsiveOnly) {
  AfirOptMockCalculator mockCalculator;
  AfirOptimizer<Lbfgs> afir(mockCalculator);
  afir.check.maxIter = 50;
  afir.optimizer.linesearch = false;
  afir.transformCoordinates = false;
  afir.check.stepMaxCoeff = 1.0e-7;
  afir.check.stepRMS = 1.0e-8;
  afir.check.gradMaxCoeff = 1.0e-7;
  afir.check.gradRMS = 1.0e-8;
  afir.check.deltaValue = 1.0e-9;
  afir.attractive = false;
  afir.phaseIn = 0;
  afir.rhsList = {0};
  afir.lhsList = {1};
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
  EXPECT_TRUE(r12 > sqrt(10));
  EXPECT_NEAR(r23, sqrt(10), 1.0e-8);
}

TEST(AfirOptimizerTests, Lbfgs_ListAndWeak) {
  AfirOptMockCalculator mockCalculator;
  AfirOptimizer<Lbfgs> afir(mockCalculator);
  afir.check.maxIter = 50;
  afir.optimizer.linesearch = false;
  afir.check.stepMaxCoeff = 1.0e-7;
  afir.check.stepRMS = 1.0e-8;
  afir.check.gradMaxCoeff = 1.0e-7;
  afir.check.gradRMS = 1.0e-8;
  afir.check.deltaValue = 1.0e-9;
  afir.phaseIn = 0;
  afir.weak = true;
  afir.attractive = false;
  afir.transformCoordinates = false;
  afir.rhsList = {0};
  afir.lhsList = {1};
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
  EXPECT_TRUE(r12 > sqrt(10));
  EXPECT_TRUE(r23 < sqrt(10));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
