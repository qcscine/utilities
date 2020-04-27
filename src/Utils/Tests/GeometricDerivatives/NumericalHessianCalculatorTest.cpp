/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/GeometricDerivatives/NumericalHessianCalculator.h"
#include "Utils/GeometryOptimization/IrcOptimizer.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

// Define a mock state
class MockState : public Core::State {};

// Define a mock calculator
class HessianMockCalculator : public CloneInterface<HessianMockCalculator, Core::Calculator> {
 public:
  HessianMockCalculator() = default;
  ~HessianMockCalculator() final = default;
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
    r_ = Results();
    auto pos1 = structure_.getPosition(0);
    auto pos2 = structure_.getPosition(1);
    double dist = (pos1 - pos2).norm();
    r_.set<Property::Energy>(cos(dist));
    GradientCollection g(structure_.size(), 3);
    g(0, 0) = -sin(dist) * (pos1[0] - pos2[0]) / dist;
    g(0, 1) = -sin(dist) * (pos1[1] - pos2[1]) / dist;
    g(0, 2) = -sin(dist) * (pos1[2] - pos2[2]) / dist;
    g(1, 0) = -g(0, 0);
    g(1, 1) = -g(0, 1);
    g(1, 2) = -g(0, 2);
    r_.set<Property::Gradients>(g);

    double charge1 = -0.5;
    double charge2 = 0.5;
    Dipole d = pos1 * charge1 + pos2 * charge2;
    r_.set<Property::Dipole>(d);

    return r_;
  };
  std::string name() const final {
    return std::string("MockCalculator");
  };
  const Settings& settings() const final {
    return settings_;
  }
  Settings& settings() final {
    return settings_;
  }
  std::unique_ptr<Utils::AtomCollection> getStructure() const final {
    return nullptr;
  }

  Utils::Results& results() final {
    return r_;
  }
  const Utils::Results& results() const final {
    return r_;
  }
  bool supportsMethodFamily(const std::string& methodFamily) const final {
    return true;
  }
  std::shared_ptr<Core::State> getState() const final {
    auto ret = std::make_shared<MockState>();
    return ret;
  }

  void loadState(std::shared_ptr<Core::State> state) final {
    auto mockstate = std::dynamic_pointer_cast<MockState>(state);
  }

 private:
  AtomCollection structure_;
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
  HessianMockCalculator mockCalculator;
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
  HessianMockCalculator mockCalculator;
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
      ASSERT_NEAR(val, results.get<Property::Hessian>()(i, j), 1.0e-5);
    }
  }
}

TEST(NumericalHessianCalculator, SeminumericalWithGradient) {
  double charge1 = -0.5;
  double charge2 = 0.5;
  HessianMockCalculator mockCalculator;
  PositionCollection pos(2, 3);
  pos(0, 0) = -M_PI / (2.0 * sqrt(3.0));
  pos(0, 1) = -M_PI / (2.0 * sqrt(3.0));
  pos(0, 2) = -M_PI / (2.0 * sqrt(3.0));
  pos(1, 0) = +M_PI / (2.0 * sqrt(3.0));
  pos(1, 1) = +M_PI / (2.0 * sqrt(3.0));
  pos(1, 2) = +M_PI / (2.0 * sqrt(3.0));
  auto elements = ElementTypeCollection{ElementType::H, ElementType::F};
  mockCalculator.setStructure(AtomCollection(elements, pos));
  NumericalHessianCalculator hessianCalc(mockCalculator);
  hessianCalc.requiredDipoleGradient(true);
  auto results = hessianCalc.calculateFromGradientDifferences(0.005);

  DipoleGradient dipoleGradientAnalytical = DipoleGradient::Zero(6, 3);
  dipoleGradientAnalytical.row(0) = Dipole{charge1, 0, 0};
  dipoleGradientAnalytical.row(1) = Dipole{0, charge1, 0};
  dipoleGradientAnalytical.row(2) = Dipole{0, 0, charge1};
  dipoleGradientAnalytical.row(3) = Dipole{charge2, 0, 0};
  dipoleGradientAnalytical.row(4) = Dipole{0, charge2, 0};
  dipoleGradientAnalytical.row(5) = Dipole{0, 0, charge2};

  DipoleGradient dipoleGradient = results.get<Utils::Property::DipoleGradient>();

  for (int dimension = 0; dimension < 6; ++dimension) {
    EXPECT_THAT(dipoleGradientAnalytical.row(dimension).x(), DoubleNear(dipoleGradient.row(dimension).x(), 1e-3));
    EXPECT_THAT(dipoleGradientAnalytical.row(dimension).y(), DoubleNear(dipoleGradient.row(dimension).y(), 1e-3));
    EXPECT_THAT(dipoleGradientAnalytical.row(dimension).z(), DoubleNear(dipoleGradient.row(dimension).z(), 1e-3));
  }
}
} // namespace Tests
} // namespace Utils
} // namespace Scine
