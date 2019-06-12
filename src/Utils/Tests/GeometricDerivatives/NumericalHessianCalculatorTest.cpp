/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/GeometricDerivatives/NumericalHessianCalculator.h"
#include "Utils/GeometryOptimization/IRCOptimizer.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

// Define a mock states handler
class HessianMockStatesHandler : public StatesHandler {
 public:
  std::shared_ptr<State> getCurrentState(StateSize /*size*/) const override {
    return nullptr;
  };
};

// Define a mock calculator
class HessianMockCalculator : public CloneInterface<HessianMockCalculator, Core::Calculator> {
 public:
  HessianMockCalculator(HessianMockStatesHandler sh) : sh_(sh){};
  ~HessianMockCalculator() = default;
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
    r_ = Results();
    auto pos1 = structure_.getPosition(0);
    auto pos2 = structure_.getPosition(1);
    double dist = (pos1 - pos2).norm();
    r_.setEnergy(cos(dist));
    GradientCollection g(structure_.size(), 3);
    g(0, 0) = -sin(dist) * (pos1[0] - pos2[0]) / dist;
    g(0, 1) = -sin(dist) * (pos1[1] - pos2[1]) / dist;
    g(0, 2) = -sin(dist) * (pos1[2] - pos2[2]) / dist;
    g(1, 0) = -g(0, 0);
    g(1, 1) = -g(0, 1);
    g(1, 2) = -g(0, 2);
    r_.setGradients(g);
    return r_;
  };
  std::string name() const override {
    return std::string("MockCalculator");
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
  std::unique_ptr<Utils::AtomCollection> getStructure() const final {
    return nullptr;
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

 private:
  AtomCollection structure_;
  HessianMockStatesHandler sh_;
  Results r_;
  Settings settings_ = Settings("dummy");
  //   Core::Calculator* cloneImpl() const override {
  //     return nullptr;
  //   }
};

/**
 * @class Scine::Utils::Tests::NumericalHessianCalculatorTest
 * @brief Comprises tests for the class Scine::Utils::NumericalHessianCalculator.
 * @test
 */
TEST(NumericalHessianCalculatorTest, Numerical) {
  HessianMockStatesHandler mockStatesHandler;
  HessianMockCalculator mockCalculator(mockStatesHandler);
  PositionCollection pos(2, 3);
  pos(0, 0) = -M_PI / (2.0 * sqrt(3.0));
  pos(0, 1) = -M_PI / (2.0 * sqrt(3.0));
  pos(0, 2) = -M_PI / (2.0 * sqrt(3.0));
  pos(1, 0) = +M_PI / (2.0 * sqrt(3.0));
  pos(1, 1) = +M_PI / (2.0 * sqrt(3.0));
  pos(1, 2) = +M_PI / (2.0 * sqrt(3.0));
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H};
  mockCalculator.setStructure(AtomCollection(elements, pos));
  NumericalHessianCalculator hessianCalc(mockCalculator);
  auto hessian = hessianCalc.calculateFromEnergyDifferences(0.005);
  for (unsigned int i = 0; i < 6; i++) {
    for (unsigned int j = 0; j < 6; j++) {
      const double val = (1.0 / 3.0) * (i < 3 ? 1.0 : -1.0) * (j < 3 ? 1.0 : -1.0);
      ASSERT_NEAR(val, hessian(i, j), 1.0e-5);
    }
  }
}

TEST(NumericalHessianCalculatorTest, SemiNumerical) {
  HessianMockStatesHandler mockStatesHandler;
  HessianMockCalculator mockCalculator(mockStatesHandler);
  PositionCollection pos(2, 3);
  pos(0, 0) = -M_PI / (2.0 * sqrt(3.0));
  pos(0, 1) = -M_PI / (2.0 * sqrt(3.0));
  pos(0, 2) = -M_PI / (2.0 * sqrt(3.0));
  pos(1, 0) = +M_PI / (2.0 * sqrt(3.0));
  pos(1, 1) = +M_PI / (2.0 * sqrt(3.0));
  pos(1, 2) = +M_PI / (2.0 * sqrt(3.0));
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H};
  mockCalculator.setStructure(AtomCollection(elements, pos));
  NumericalHessianCalculator hessianCalc(mockCalculator);
  auto results = hessianCalc.calculateFromGradientDifferences(0.005);
  for (unsigned int i = 0; i < 6; i++) {
    for (unsigned int j = 0; j < 6; j++) {
      const double val = (1.0 / 3.0) * (i < 3 ? 1.0 : -1.0) * (j < 3 ? 1.0 : -1.0);
      ASSERT_NEAR(val, results.getHessian()(i, j), 1.0e-5);
    }
  }
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
