/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/GeometryOptimization/NtOptimizer.h"
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
class NtOptMockCalculator : public CloneInterface<NtOptMockCalculator, Core::Calculator> {
 public:
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
    auto v12 = p1 - p2;
    auto r12 = v12.norm();

    // Energy
    // f(r12) = (r12-8)(r12-8)(r12-4)(r12-12)
    double e = (r12 - 8.0) * (r12 - 8.0) * (r12 - 4.0) * (r12 - 12.0);
    // Gradient
    GradientCollection g(structure_.size(), 3);
    g(0, 0) = (4.0 * (pow(r12, 3.0) - 24.0 * pow(r12, 2.0) + 184.0 * r12 - 448.0) / r12) * v12[0];
    g(0, 1) = (4.0 * (pow(r12, 3.0) - 24.0 * pow(r12, 2.0) + 184.0 * r12 - 448.0) / r12) * v12[1];
    g(0, 2) = (4.0 * (pow(r12, 3.0) - 24.0 * pow(r12, 2.0) + 184.0 * r12 - 448.0) / r12) * v12[2];
    g(1, 0) = -(4.0 * (pow(r12, 3.0) - 24.0 * pow(r12, 2.0) + 184.0 * r12 - 448.0) / r12) * v12[0];
    g(1, 1) = -(4.0 * (pow(r12, 3.0) - 24.0 * pow(r12, 2.0) + 184.0 * r12 - 448.0) / r12) * v12[1];
    g(1, 2) = -(4.0 * (pow(r12, 3.0) - 24.0 * pow(r12, 2.0) + 184.0 * r12 - 448.0) / r12) * v12[2];

    r_ = Results();
    r_.set<Property::SuccessfulCalculation>(true);
    r_.set<Property::Energy>(e);
    r_.set<Property::Gradients>(g);
    return r_;
  };
  std::string name() const final {
    return std::string("NtOptMockCalculator");
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

// Define a test logger
namespace {
Core::Log logger = Core::Log::silent();
}

/**
 * @class Scine::Utils::Tests::NtOptimizerTests
 * @brief Comprises tests for the class Scine::Utils::NtOptimizer.
 * @test
 */
TEST(NtOptimizerTests, DefaultAttraction) {
  NtOptMockCalculator mockCalculator;
  NtOptimizer nt(mockCalculator);
  nt.lhsList = {0};
  nt.rhsList = {1};
  nt.totalForceNorm = 1.0;
  nt.check.maxIter = 600;
  nt.optimizer.factor = 0.02;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(2, 3);
  positions(0, 1) = 11.0;
  AtomCollection atoms(elements, positions);
  auto nIter = nt.optimize(atoms, logger);
  // Check results
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto r12 = (p1 - p2).norm();
  EXPECT_TRUE(nt.check.maxIter > nIter);
  EXPECT_NEAR(r12, 8.0, 1.0e-2);
}

TEST(NtOptimizerTests, DefaultRepulsion) {
  NtOptMockCalculator mockCalculator;
  NtOptimizer nt(mockCalculator);
  nt.lhsList = {0};
  nt.rhsList = {1};
  nt.attractive = false;
  nt.check.repulsiveStop = 10.0;
  nt.totalForceNorm = 1.0;
  nt.check.maxIter = 400;
  nt.optimizer.factor = 0.02;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(2, 3);
  positions(0, 1) = 5.0;
  AtomCollection atoms(elements, positions);
  auto nIter = nt.optimize(atoms, logger);
  // Check results
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto r12 = (p1 - p2).norm();
  EXPECT_TRUE(nt.check.maxIter > nIter);
  EXPECT_NEAR(r12, 8.0, 1.0e-2);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
