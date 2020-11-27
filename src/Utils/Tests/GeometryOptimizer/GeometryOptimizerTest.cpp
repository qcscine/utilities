/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/GeometryOptimization/GeometryOptimizer.h"
#include "Utils/GeometryOptimization/IrcOptimizer.h"
#include "Utils/Optimizer/GradientBased/Lbfgs.h"
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

// Define a mock calculator
class GeoOptMockCalculator : public Core::Calculator {
 public:
  GeoOptMockCalculator() = default;
  ~GeoOptMockCalculator() = default;
  void setStructure(const AtomCollection& structure) final {
    structure_ = structure;
  }
  void modifyPositions(PositionCollection newPositions) final {
    structure_.setPositions(newPositions);
  }
  const PositionCollection& getPositions() const final {
    return structure_.getPositions();
  }
  void setRequiredProperties(const PropertyList& requiredProperties) final{};
  PropertyList getRequiredProperties() const final {
    return Utils::PropertyList{};
  }
  PropertyList possibleProperties() const final {
    return Utils::Property::Energy | Utils::Property::Gradients;
  }
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

    // Hessian
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(9, 9);
    // diagonal
    h(0, 0) = 8 * v12[0] * v12[0] - (4 * pow(r12, 2) - 40);
    h(1, 1) = 8 * v12[1] * v12[1] - (4 * pow(r12, 2) - 40);
    h(2, 2) = 8 * v12[2] * v12[2] - (4 * pow(r12, 2) - 40);
    h(6, 6) = 8 * v32[0] * v32[0] - (4 * pow(r23, 2) - 40);
    h(7, 7) = 8 * v32[1] * v32[1] - (4 * pow(r23, 2) - 40);
    h(8, 8) = 8 * v32[2] * v32[2] - (4 * pow(r23, 2) - 40);
    h(3, 3) = h(0, 0) + h(6, 6);
    h(4, 4) = h(1, 1) + h(7, 7);
    h(5, 5) = h(2, 2) + h(8, 8);
    // rest of h-1-1
    h(1, 0) = 8 * v12[1] * v12[0];
    h(2, 0) = 8 * v12[2] * v12[0];
    h(2, 1) = 8 * v12[2] * v12[1];
    h(0, 1) = h(1, 0);
    h(0, 2) = h(2, 0);
    h(1, 2) = h(2, 1);
    // rest of h-3-3
    h(7, 6) = 8 * v32[1] * v32[0];
    h(8, 6) = 8 * v32[2] * v32[0];
    h(8, 7) = 8 * v32[2] * v32[1];
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
    r_.set<Property::SuccessfulCalculation>(true);
    r_.set<Property::Energy>(e);
    r_.set<Property::Gradients>(g);
    r_.set<Property::Hessian>(h);
    return r_;
  };
  std::string name() const final {
    return std::string("GeoOptMockCalculator");
  };
  std::shared_ptr<Core::State> getState() const final {
    return nullptr;
  }
  void loadState(std::shared_ptr<Core::State> state) final {
  }
  const Settings& settings() const final {
    return settings_;
  }
  Settings& settings() final {
    return settings_;
  }
  virtual Utils::Results& results() final {
    return r_;
  }
  virtual const Utils::Results& results() const final {
    return r_;
  }
  virtual std::unique_ptr<Utils::AtomCollection> getStructure() const final {
    return nullptr;
  }
  bool supportsMethodFamily(const std::string& methodFamily) const final {
    return true;
  }

 private:
  AtomCollection structure_;
  Results r_;
  Settings settings_ = Settings("dummy");
  Core::Calculator* cloneImpl() const final {
    return nullptr;
  }
};

// Define a test logger
namespace {
Core::Log logger = Core::Log::silent();
}

/**
 * @class Scine::Utils::Tests::GeometryOptimizerTests
 * @brief Comprises tests for the class Scine::Utils::GeometryOptimizer.
 * @test
 */
TEST(GeometryOptimizerTests, SteepestDescent) {
  GeoOptMockCalculator mockCalculator;
  GeometryOptimizer<SteepestDescent> geo(mockCalculator);
  geo.transformCoordinates = false;
  geo.check.maxIter = 50;
  geo.optimizer.factor = 0.01;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-9;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-9;
  geo.check.deltaValue = 1.0e-12;
  geo.check.requirement = 4;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -M_PI;
  positions(2, 0) = M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = geo.optimize(atoms, logger);

  // Check results
  EXPECT_TRUE(nIter < geo.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, sqrt(10), 1.0e-8);
  EXPECT_NEAR(r23, sqrt(10), 1.0e-8);
}

TEST(GeometryOptimizerTests, Bfgs) {
  GeoOptMockCalculator mockCalculator;
  GeometryOptimizer<Bfgs> geo(mockCalculator);
  geo.transformCoordinates = false;
  geo.optimizer.useTrustRadius = true;
  geo.check.maxIter = 20;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-8;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-8;
  geo.check.deltaValue = 1.0e-12;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = 0.9 * M_PI;
  positions(2, 0) = 0.9 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = geo.optimize(atoms, logger);

  // Check results
  EXPECT_TRUE(nIter < geo.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, sqrt(10), 1.0e-8);
  EXPECT_NEAR(r23, sqrt(10), 1.0e-8);
}

TEST(GeometryOptimizerTests, Bfgs_NoDiis) {
  GeoOptMockCalculator mockCalculator;
  GeometryOptimizer<Bfgs> geo(mockCalculator);
  geo.transformCoordinates = false;
  geo.optimizer.useGdiis = false;
  geo.optimizer.useTrustRadius = true;
  geo.check.maxIter = 30;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-8;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-8;
  geo.check.deltaValue = 1.0e-12;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = 0.75 * M_PI;
  positions(2, 0) = 0.64 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = geo.optimize(atoms, logger);

  // Check results
  EXPECT_TRUE(nIter < geo.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, sqrt(10), 1.0e-8);
  EXPECT_NEAR(r23, sqrt(10), 1.0e-8);
}

TEST(GeometryOptimizerTests, Lbfgs) {
  GeoOptMockCalculator mockCalculator;
  GeometryOptimizer<Lbfgs> geo(mockCalculator);
  geo.transformCoordinates = false;
  geo.optimizer.maxBacktracking = 20;
  geo.check.maxIter = 300;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-8;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-8;
  geo.check.deltaValue = 1.0e-12;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = 0.75 * M_PI;
  positions(2, 0) = 0.64 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = geo.optimize(atoms, logger);

  // Check results
  EXPECT_TRUE(nIter < geo.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, sqrt(10), 1.0e-8);
  EXPECT_NEAR(r23, sqrt(10), 1.0e-8);
}

TEST(GeometryOptimizerTests, LbfgsNoLineSearch) {
  GeoOptMockCalculator mockCalculator;
  GeometryOptimizer<Lbfgs> geo(mockCalculator);
  geo.transformCoordinates = false;
  geo.optimizer.linesearch = false;
  geo.check.maxIter = 20;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-9;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-9;
  geo.check.deltaValue = 1.0e-12;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = 0.75 * M_PI;
  positions(2, 0) = 0.64 * M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = geo.optimize(atoms, logger);

  // Check results
  EXPECT_TRUE(nIter < geo.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, sqrt(10), 1.0e-8);
  EXPECT_NEAR(r23, sqrt(10), 1.0e-8);
}

TEST(GeometryOptimizerTests, NewtonRaphson) {
  GeoOptMockCalculator mockCalculator;
  GeometryOptimizer<NewtonRaphson> geo(mockCalculator);
  geo.check.maxIter = 50;
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
  auto nIter = geo.optimize(atoms, logger);

  // Check results
  EXPECT_TRUE(nIter < geo.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, sqrt(10), 1.0e-8);
  EXPECT_NEAR(r23, sqrt(10), 1.0e-8);
}

TEST(GeometryOptimizerTests, TestObserver) {
  GeoOptMockCalculator mockCalculator;
  GeometryOptimizer<SteepestDescent> geo(mockCalculator);
  geo.optimizer.factor = 0.01;
  geo.check.maxIter = 63;
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
  auto nIter = geo.optimize(atoms, logger);
  // Check
  EXPECT_EQ(nIter, counter);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
