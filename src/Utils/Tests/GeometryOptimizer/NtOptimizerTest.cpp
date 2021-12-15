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
#include <Utils/CalculatorBasics/TestCalculator.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Math/QuaternionFit.h>
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
    r_ = Results{};
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
    return Utils::PropertyList{};
  };
  PropertyList possibleProperties() const final {
    return Utils::Property::Energy | Utils::Property::Gradients;
  };
  const Results& calculate(std::string /*dummy*/ = "") final {
    auto p1 = structure_.getPosition(0);
    auto p2 = structure_.getPosition(1);
    auto v12 = p1 - p2;
    auto r12 = v12.norm();

    // Energy
    // f(r12) = (r12-8)(r12-8)(r12-4)(r12-12)
    double e = (r12 - 8.0) * (r12 - 8.0) * (r12 - 4.0) * (r12 - 12.0);
    // Gradient
    GradientCollection g(structure_.size(), 3);
    g(0, 0) = (4.0 * (std::pow(r12, 3.0) - 24.0 * std::pow(r12, 2.0) + 184.0 * r12 - 448.0) / r12) * v12[0];
    g(0, 1) = (4.0 * (std::pow(r12, 3.0) - 24.0 * std::pow(r12, 2.0) + 184.0 * r12 - 448.0) / r12) * v12[1];
    g(0, 2) = (4.0 * (std::pow(r12, 3.0) - 24.0 * std::pow(r12, 2.0) + 184.0 * r12 - 448.0) / r12) * v12[2];
    g(1, 0) = -(4.0 * (std::pow(r12, 3.0) - 24.0 * std::pow(r12, 2.0) + 184.0 * r12 - 448.0) / r12) * v12[0];
    g(1, 1) = -(4.0 * (std::pow(r12, 3.0) - 24.0 * std::pow(r12, 2.0) + 184.0 * r12 - 448.0) / r12) * v12[1];
    g(1, 2) = -(4.0 * (std::pow(r12, 3.0) - 24.0 * std::pow(r12, 2.0) + 184.0 * r12 - 448.0) / r12) * v12[2];

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
  void loadState(std::shared_ptr<Core::State> /* state */) final {
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
  bool supportsMethodFamily(const std::string& /*methodFamily*/) const final {
    return true;
  }

 private:
  AtomCollection structure_;
  Results r_;
  Settings settings_ = Settings("dummy");
};

/**
 * @class Scine::Utils::Tests::NtOptimizerTests
 * @brief Comprises tests for the class Scine::Utils::NtOptimizer.
 * @test
 */
TEST(NtOptimizerTests, DefaultAttraction) {
  Core::Log logger = Core::Log::silent();
  NtOptMockCalculator mockCalculator;
  NtOptimizer nt(mockCalculator);
  nt.lhsList = {0};
  nt.rhsList = {1};
  nt.totalForceNorm = 1.0;
  nt.check.maxIter = 500;
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
  EXPECT_TRUE(nt.check.maxIter > static_cast<unsigned>(nIter));
  EXPECT_EQ(nIter, 497);
  EXPECT_NEAR(r12, 8.0, 1.0e-8);
}

TEST(NtOptimizerTests, SN2Attraction) {
  Core::Log logger = Core::Log::silent();
  TestCalculator testCalc;
  NtOptimizer nt(testCalc);
  nt.coordinateSystem = CoordinateSystem::Cartesian;
  nt.lhsList = {0};
  nt.rhsList = {1};
  nt.numberOfMicroCycles = 10;
  auto elements = ElementTypeCollection{ElementType::Cl, ElementType::C, ElementType::H,
                                        ElementType::H,  ElementType::H, ElementType::Br};
  PositionCollection startPositions = PositionCollection::Zero(6, 3);
  // clang-format off
  startPositions <<  3.7376961460e+00,  2.1020866350e-04,  4.5337439168e-02,
                    -3.8767703481e+00, -2.4803422157e-05, -1.2049608882e-01,
                    -2.3620148614e+00,  1.3238308540e+00,  1.0376490681e-01,
                    -2.3809041075e+00, -8.2773666259e-01,  9.6331578315e-01,
                    -2.3309449521e+00, -4.9652606314e-01, -1.3293307598e+00,
                    -7.4798903722e+00,  2.6536371103e-04, -1.9897114399e-01;
  // clang-format on
  AtomCollection atoms(elements, startPositions);
  auto nIter = nt.optimize(atoms, logger);
  ASSERT_EQ(nIter, 48);
  // Check results
  PositionCollection reference = PositionCollection::Zero(6, 3);
  // clang-format off
  reference <<  1.9881110238e+00,  1.5620965399e-04,  7.2336700390e-03,
               -2.1271852259e+00,  2.9195587350e-05, -8.2392319691e-02,
               -1.7301242740e+00,  1.9682197934e+00,  2.1134360238e-01,
               -1.7588564787e+00, -1.2334301806e+00,  1.4862704381e+00,
               -1.6845481419e+00, -7.4219949889e-01, -1.9184062691e+00,
               -5.9073039827e+00,  3.6871285900e-03, -1.6475583032e-01;
  // clang-format on
  PositionCollection positions = atoms.getPositions();
  for (unsigned int i = 0; i < reference.array().size(); i++) {
    EXPECT_NEAR(reference.data()[i], positions.data()[i], 1.0e-8);
  }
  // do again with one side only
  atoms.setPositions(startPositions);
  nt.movableSide = "rhs";
  nIter = nt.optimize(atoms, logger);
  positions = atoms.getPositions();
  auto fit = Utils::QuaternionFit(reference, positions);
  EXPECT_LT(fit.getRMSD(), 0.05);
  ASSERT_EQ(nIter, 48);
}

TEST(NtOptimizerTests, AttractionSingleStep) {
  Core::Log logger = Core::Log::silent();
  NtOptMockCalculator mockCalculator;
  NtOptimizer nt(mockCalculator);
  nt.lhsList = {0};
  nt.rhsList = {1};
  nt.totalForceNorm = 1.0;
  nt.check.maxIter = 10;
  nt.optimizer.factor = 0.02;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(2, 3);
  // Enforce a scan that stops after the initial evaluation.
  positions(0, 1) = 0.1;
  AtomCollection atoms(elements, positions);
  ASSERT_THROW(nt.optimize(atoms, logger), std::runtime_error);
}

TEST(NtOptimizerTests, DefaultRepulsion) {
  Core::Log logger = Core::Log::silent();
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
  EXPECT_TRUE(nt.check.maxIter > static_cast<unsigned>(nIter));
  EXPECT_EQ(nIter, 356);
  EXPECT_NEAR(r12, 8.0, 1.0e-8);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
