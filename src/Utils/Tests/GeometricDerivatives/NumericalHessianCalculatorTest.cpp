/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/GeometricDerivatives/NumericalHessianCalculator.h"
#include "Utils/GeometryOptimization/IrcOptimizer.h"
#include "Utils/MSVCCompatibility.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>
#include <Utils/Math/AutomaticDifferentiation/VectorDerivatives3D.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <gmock/gmock.h>
#include <omp.h>
#include <vector>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

// Define a mock state
class MockState final : public Core::State {};

// Define a mock calculator
class HessianMockCalculator : public CloneInterface<HessianMockCalculator, Core::Calculator> {
 public:
  void setStructure(const AtomCollection& structure) final {
    atomicCharges_ = Eigen::VectorXd::Random(structure.size());
    atomicCharges_(0) = -atomicCharges_.sum() + atomicCharges_(0);
    r_ = Results();
    r_.set<Property::AtomicCharges>(std::vector<double>(atomicCharges_.data(), atomicCharges_.data() + atomicCharges_.size()));
    structure_ = structure;
  };
  void modifyPositions(PositionCollection newPositions) final {
    structure_.setPositions(newPositions);
    r_ = Results{};
  };
  const PositionCollection& getPositions() const final {
    return structure_.getPositions();
  };
  void setRequiredProperties(const PropertyList& /* requiredProperties */) final{};
  PropertyList getRequiredProperties() const final {
    return PropertyList{};
  };
  PropertyList possibleProperties() const final {
    return Property::Energy | Property::Gradients | Property::Dipole | Property::DipoleGradient | Property::AtomicCharges;
  };
  const Results& calculate(std::string /*dummy*/) final {
    double energy = 0.0;
    Dipole d = Dipole::Zero();
    DipoleGradient dipGrad = DipoleGradient::Zero(getPositions().size(), 3);
    auto dip = AutomaticDifferentiation::VectorDerivatives3D::spatialVectorHessian3D(d);
    for (int atom_i = 0; atom_i < getPositions().rows(); ++atom_i) {
      auto atomicDipContrib =
          AutomaticDifferentiation::VectorDerivatives3D::spatialVectorHessian3D(getPositions().row(atom_i)) *
          atomicCharges_[atom_i];
      dip = dip + atomicDipContrib;
      dipGrad.block<3, 1>(atom_i * 3, 0) += Gradient(atomicDipContrib.x().deriv());
      dipGrad.block<3, 1>(atom_i * 3, 1) += Gradient(atomicDipContrib.y().deriv());
      dipGrad.block<3, 1>(atom_i * 3, 2) += Gradient(atomicDipContrib.z().deriv());
    }

    GradientCollection g = GradientCollection::Zero(structure_.size(), 3);
    for (int atom_i = 0; atom_i < getPositions().rows(); ++atom_i) {
      for (int atom_j = atom_i; atom_j < getPositions().rows(); ++atom_j) {
        Position distanceVector = getPositions().row(atom_j) - getPositions().row(atom_i);
        auto distance = AutomaticDifferentiation::variableWithUnitDerivative<DerivativeOrder::One>(distanceVector.norm());
        auto energyWithDeriv = cos(distance);
        energy += energyWithDeriv.value();
        auto gradient = AutomaticDifferentiation::get3Dfrom1D<DerivativeOrder::One>(energyWithDeriv, distanceVector);
        AutomaticDifferentiation::addDerivativeToContainer<Derivative::First>(g, atom_i, atom_j,
                                                                              Gradient(gradient.derivatives()));
      }
    }

#pragma omp critical(writeData)
    {
      r_.set<Property::Energy>(energy);
      r_.set<Property::Dipole>({dip.x().value(), dip.y().value(), dip.z().value()});
      r_.set<Property::Gradients>(g);
      r_.set<Property::DipoleGradient>(dipGrad);
      r_.set<Property::SuccessfulCalculation>(true);
    }

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
  bool supportsMethodFamily(const std::string& /* methodFamily */) const final {
    return true;
  }
  std::shared_ptr<Core::State> getState() const final {
    auto ret = std::make_shared<MockState>();
    return ret;
  }

  void loadState(std::shared_ptr<Core::State> state) final {
    auto mockstate = std::dynamic_pointer_cast<MockState>(state);
  }
  bool allowsPythonGILRelease() const override {
    return true;
  };

 private:
  AtomCollection structure_;
  Eigen::VectorXd atomicCharges_;
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
  auto results2 = hessianCalc.calculate(0.005);
  for (unsigned int i = 0; i < 6; i++) {
    for (unsigned int j = 0; j < 6; j++) {
      const double val = (1.0 / 3.0) * (i < 3 ? 1.0 : -1.0) * (j < 3 ? 1.0 : -1.0);
      EXPECT_NEAR(val, results.get<Property::Hessian>()(i, j), 1.0e-5);
      EXPECT_NEAR(results.get<Property::Hessian>()(i, j), results2.get<Property::Hessian>()(i, j), 1.0e-5);
    }
  }
}

TEST(NumericalHessianCalculator, SeminumericalWithGradient) {
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
  auto charges = mockCalculator.results().get<Property::AtomicCharges>();

  mockCalculator.calculate("");
  DipoleGradient dipoleGradientAnalytical = mockCalculator.results().get<Property::DipoleGradient>();

  DipoleGradient dipoleGradient = results.get<Utils::Property::DipoleGradient>();

  for (int dimension = 0; dimension < 6; ++dimension) {
    EXPECT_THAT(dipoleGradientAnalytical.row(dimension).x(), DoubleNear(dipoleGradient.row(dimension).x(), 1e-3));
    EXPECT_THAT(dipoleGradientAnalytical.row(dimension).y(), DoubleNear(dipoleGradient.row(dimension).y(), 1e-3));
    EXPECT_THAT(dipoleGradientAnalytical.row(dimension).z(), DoubleNear(dipoleGradient.row(dimension).z(), 1e-3));
  }
}

TEST(NumericalHessianCalculator, SeminumericalWithGradientAndSelectedAtoms) {
  HessianMockCalculator mockCalculator;
  PositionCollection pos(3, 3);
  pos(0, 0) = -M_PI / (2.0 * sqrt(3.0));
  pos(0, 1) = -M_PI / (2.0 * sqrt(3.0));
  pos(0, 2) = -M_PI / (2.0 * sqrt(3.0));
  pos(1, 0) = +4. * M_PI / (2.0 * sqrt(3.0));
  pos(1, 1) = +4. * M_PI / (2.0 * sqrt(3.0));
  pos(1, 2) = +4. * M_PI / (2.0 * sqrt(3.0));
  pos(2, 0) = +5. * M_PI / (2.0 * sqrt(3.0));
  pos(2, 1) = +5. * M_PI / (2.0 * sqrt(3.0));
  pos(2, 2) = +5. * M_PI / (2.0 * sqrt(3.0));
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::F};
  mockCalculator.setStructure(AtomCollection(elements, pos));
  NumericalHessianCalculator hessianCalc(mockCalculator);
  hessianCalc.requiredDipoleGradient(true);
  auto results = hessianCalc.calculateFromGradientDifferences(0.005);
  const HessianMatrix& hessian1 = results.get<Property::Hessian>();

  mockCalculator.calculate("");
  DipoleGradient dipoleGradientAnalytical = mockCalculator.results().get<Property::DipoleGradient>();

  auto indices = std::vector<int>(2);
  indices[0] = 0;
  indices[1] = 2;
  auto results2 = hessianCalc.calculateFromGradientDifferences(0.005, indices);
  auto results3 = hessianCalc.calculate(indices, 0.005);
  const HessianMatrix& hessian2 = results2.get<Property::Hessian>();
  const HessianMatrix& hessian3 = results3.get<Property::Hessian>();
  DipoleGradient dipoleGradient = results2.get<Utils::Property::DipoleGradient>();

  for (int row = 0; row < hessian1.rows(); ++row) {
    for (int col = 0; col < hessian1.cols(); ++col) {
      if (hessian2(row, col) > 1.0e-13) {
        EXPECT_NEAR(hessian2(row, col), hessian1(row, col), 1e-13);
        EXPECT_NEAR(hessian3(row, col), hessian1(row, col), 1e-13);
      }
    }
  }

  for (int index : indices) {
    EXPECT_NEAR(dipoleGradient.row(index).x(), dipoleGradientAnalytical.row(index).x(), 1.0e-13);
    EXPECT_NEAR(dipoleGradient.row(index).y(), dipoleGradientAnalytical.row(index).y(), 1.0e-13);
    EXPECT_NEAR(dipoleGradient.row(index).z(), dipoleGradientAnalytical.row(index).z(), 1.0e-13);
  }
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
