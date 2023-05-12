/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/GeometryOptimization/GeometryOptimization.h"
#include "Utils/GeometryOptimization/IrcOptimizer.h"
#include "Utils/GeometryOptimization/UnitCellGeometryOptimizer.h"
#include "Utils/LennardJonesCalculator/LennardJonesCalculatorSettings.h"
#include "Utils/Optimizer/GradientBased/Lbfgs.h"
#include "Utils/Optimizer/GradientBased/SteepestDescent.h"
#include "Utils/Optimizer/HessianBased/EigenvectorFollowing.h"
#include "Utils/Optimizer/HessianBased/NewtonRaphson.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/CalculatorBasics/StatesHandler.h>
#include <Utils/CalculatorBasics/TestCalculator.h>
#include <Utils/Constants.h>
#include <Utils/LennardJonesCalculator/LennardJonesCalculator.h>
#include <Utils/Settings.h>
#include <Utils/Typenames.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

// Define mock calculator settings
class UnitCellMockOptSettings : public Settings {
 public:
  UnitCellMockOptSettings() : Settings("UnitCellGeoOptMockCalculatorSettings") {
    UniversalSettings::DoubleDescriptor convergence_threshold("Energy convergence limit.");
    convergence_threshold.setDefaultValue(1e-12);
    this->_fields.push_back(SettingsNames::selfConsistenceCriterion, convergence_threshold);

    UniversalSettings::StringDescriptor pbc("PBC");
    pbc.setDefaultValue("15.0,15.0,15.0,120.0,90.0,90.0");
    this->_fields.push_back(SettingsNames::periodicBoundaries, pbc);

    this->resetToDefaults();
  };
  ~UnitCellMockOptSettings() override = default;
};

// Define a mock calculator
class UnitCellGeoOptMockCalculator : public Core::Calculator {
 public:
  UnitCellGeoOptMockCalculator() {
    settings_ = std::make_unique<UnitCellMockOptSettings>();
  };
  ~UnitCellGeoOptMockCalculator() override = default;
  void setStructure(const AtomCollection& structure) final {
    structure_ = structure;
  }
  void modifyPositions(PositionCollection newPositions) final {
    structure_.setPositions(newPositions);
  }
  const PositionCollection& getPositions() const final {
    return structure_.getPositions();
  }
  void setRequiredProperties(const PropertyList& /* requiredProperties */) final{};
  PropertyList getRequiredProperties() const final {
    return Utils::PropertyList{};
  }
  PropertyList possibleProperties() const final {
    return Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::StressTensor;
  }
  const Results& calculate(std::string /*dummy*/ = "") final {
    const auto p1 = structure_.getPosition(0);
    const auto p2 = structure_.getPosition(1);
    const auto p3 = structure_.getPosition(2);
    const auto v12 = p1 - p2;
    const auto v32 = p3 - p2;
    const auto r12 = v12.norm();
    const auto r23 = v32.norm();

    // Energy
    // f(r12, r23) = (r12-4)*(r12-2)*(r12+2)*(r12+4) + (r23-4)*(r23-2)*(r23+2)*(r23+4);
    const double e = std::pow(r12, 4) - 20 * std::pow(r12, 2) + std::pow(r23, 4) - 20 * std::pow(r23, 2) + 128;
    // Gradient
    GradientCollection g(structure_.size(), 3);
    g(0, 0) = (4.0 * std::pow(r12, 2.0) - 40.0) * v12[0];
    g(0, 1) = (4.0 * std::pow(r12, 2.0) - 40.0) * v12[1];
    g(0, 2) = (4.0 * std::pow(r12, 2.0) - 40.0) * v12[2];
    g(1, 0) = -(4.0 * std::pow(r12, 2.0) - 40.0) * v12[0] - (4.0 * std::pow(r23, 2.0) - 40.0) * v32[0];
    g(1, 1) = -(4.0 * std::pow(r12, 2.0) - 40.0) * v12[1] - (4.0 * std::pow(r23, 2.0) - 40.0) * v32[1];
    g(1, 2) = -(4.0 * std::pow(r12, 2.0) - 40.0) * v12[2] - (4.0 * std::pow(r23, 2.0) - 40.0) * v32[2];
    g(2, 0) = (4.0 * std::pow(r23, 2.0) - 40.0) * v32[0];
    g(2, 1) = (4.0 * std::pow(r23, 2.0) - 40.0) * v32[1];
    g(2, 2) = (4.0 * std::pow(r23, 2.0) - 40.0) * v32[2];

    // Hessian
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(9, 9);
    // diagonal
    h(0, 0) = 8 * v12[0] * v12[0] - (4 * std::pow(r12, 2) - 40);
    h(1, 1) = 8 * v12[1] * v12[1] - (4 * std::pow(r12, 2) - 40);
    h(2, 2) = 8 * v12[2] * v12[2] - (4 * std::pow(r12, 2) - 40);
    h(6, 6) = 8 * v32[0] * v32[0] - (4 * std::pow(r23, 2) - 40);
    h(7, 7) = 8 * v32[1] * v32[1] - (4 * std::pow(r23, 2) - 40);
    h(8, 8) = 8 * v32[2] * v32[2] - (4 * std::pow(r23, 2) - 40);
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

    // stress tensor
    auto pbc = Utils::PeriodicBoundaries(settings_->getString(SettingsNames::periodicBoundaries));
    Eigen::Matrix3d s = this->getIdealPbc().getCellMatrix() - pbc.getCellMatrix();

    r_ = Results();
    r_.set<Property::SuccessfulCalculation>(true);
    r_.set<Property::Energy>(e);
    r_.set<Property::Gradients>(g);
    r_.set<Property::Hessian>(h);
    r_.set<Property::StressTensor>(s);
    return r_;
  };
  std::string name() const final {
    return std::string("UnitCellGeoOptMockCalculator");
  };
  std::shared_ptr<Core::State> getState() const final {
    return nullptr;
  }
  void loadState(std::shared_ptr<Core::State> /* state */) final {
  }
  const Settings& settings() const final {
    return *settings_;
  }
  Settings& settings() final {
    return *settings_;
  }
  Utils::Results& results() final {
    return r_;
  }
  const Utils::Results& results() const final {
    return r_;
  }
  std::unique_ptr<Utils::AtomCollection> getStructure() const final {
    return nullptr;
  }
  bool supportsMethodFamily(const std::string& /*methodFamily*/) const final {
    return true;
  }

  inline static PeriodicBoundaries getIdealPbc() {
    return PeriodicBoundaries(10.0);
  };

 private:
  AtomCollection structure_;
  Results r_;
  std::unique_ptr<Settings> settings_;
  std::shared_ptr<Core::Calculator> cloneImpl() const final {
    return nullptr;
  }
};

/**
 * @class Scine::Utils::Tests::UnitCellGeometryOptimizerTests
 * @brief Comprises tests for the class Scine::Utils::UnitCellGeometryOptimizer.
 * @test
 */
TEST(UnitCellGeometryOptimizerTests, SteepestDescent) {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<UnitCellGeoOptMockCalculator>();
  UnitCellGeometryOptimizer<SteepestDescent, SteepestDescent> geo(*mockCalculator);
  geo.coordinateSystem = CoordinateSystem::Cartesian;
  geo.check.maxIter = 35;
  geo.cellOptimizer.factor = 0.01;
  geo.geoOptimizer.optimizer.factor = 0.01;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-9;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-9;
  geo.check.deltaValue = 1.0e-12;
  geo.check.requirement = 4;
  ASSERT_TRUE(GeometryOptimization::settingsMakeSense(geo));
  const ElementTypeCollection elements{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -M_PI;
  positions(2, 0) = M_PI;
  AtomCollection atoms(elements, positions);
  auto pbc = PeriodicBoundaries(mockCalculator->settings().getString(SettingsNames::periodicBoundaries));
  auto system = PeriodicSystem(pbc, atoms);
  auto nIter = geo.optimize(system, logger);

  // Check results
  EXPECT_TRUE(static_cast<unsigned>(nIter) < geo.check.maxIter);
  auto p1 = system.atoms.getPosition(0);
  auto p2 = system.atoms.getPosition(1);
  auto p3 = system.atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, std::sqrt(10), 1.0e-8);
  EXPECT_NEAR(r23, std::sqrt(10), 1.0e-8);
  EXPECT_TRUE(system.pbc.getCellMatrix().isApprox(UnitCellGeoOptMockCalculator::getIdealPbc().getCellMatrix(), 1e-7));
}

TEST(UnitCellGeometryOptimizerTests, MixedOptimizer) {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<UnitCellGeoOptMockCalculator>();
  UnitCellGeometryOptimizer<SteepestDescent, Bfgs> geo(*mockCalculator);
  geo.coordinateSystem = CoordinateSystem::Cartesian;
  geo.check.maxIter = 300;
  geo.cellOptimizer.useTrustRadius = true;
  geo.geoOptimizer.optimizer.factor = 0.001;
  geo.geoOptimizer.optimizer.useTrustRadius = true;
  geo.geoOptimizer.optimizer.trustRadius = 0.005;
  geo.geoOptimizer.optimizer.dynamicMultiplier = 1.1;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-9;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-9;
  geo.check.deltaValue = 1.0e-12;
  geo.check.requirement = 4;
  ASSERT_TRUE(GeometryOptimization::settingsMakeSense(geo));
  const ElementTypeCollection elements{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -M_PI;
  positions(2, 0) = M_PI;
  AtomCollection atoms(elements, positions);
  auto nIter = geo.optimize(atoms, logger);

  // Check results
  EXPECT_TRUE(static_cast<unsigned>(nIter) < geo.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, std::sqrt(10), 1.0e-8);
  EXPECT_NEAR(r23, std::sqrt(10), 1.0e-8);
  EXPECT_TRUE(PeriodicBoundaries(mockCalculator->settings().getString(SettingsNames::periodicBoundaries))
                  .getCellMatrix()
                  .isApprox(UnitCellGeoOptMockCalculator::getIdealPbc().getCellMatrix(), 1e-7));
}

TEST(UnitCellGeometryOptimizerTests, Bfgs) {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<UnitCellGeoOptMockCalculator>();
  UnitCellGeometryOptimizer<Bfgs, Bfgs> geo(*mockCalculator);
  geo.coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
  geo.check.maxIter = 50;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-8;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-8;
  geo.check.deltaValue = 1.0e-12;
  ASSERT_TRUE(GeometryOptimization::settingsMakeSense(geo));
  const ElementTypeCollection elements{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = 0.9 * M_PI;
  positions(2, 0) = 0.9 * M_PI;
  AtomCollection atoms(elements, positions);
  auto pbc = PeriodicBoundaries(mockCalculator->settings().getString(SettingsNames::periodicBoundaries));
  auto nIter = geo.optimize(atoms, pbc, logger);

  // Check results
  EXPECT_TRUE(static_cast<unsigned>(nIter) < geo.check.maxIter);
  auto p1 = atoms.getPosition(0);
  auto p2 = atoms.getPosition(1);
  auto p3 = atoms.getPosition(2);
  auto r12 = (p1 - p2).norm();
  auto r23 = (p3 - p2).norm();
  EXPECT_NEAR(r12, std::sqrt(10), 1.0e-8);
  EXPECT_NEAR(r23, std::sqrt(10), 1.0e-8);
  EXPECT_TRUE(pbc.getCellMatrix().isApprox(UnitCellGeoOptMockCalculator::getIdealPbc().getCellMatrix(), 1e-7));
}

TEST(UnitCellGeometryOptimizerTests, TestObserver) {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<UnitCellGeoOptMockCalculator>();
  UnitCellGeometryOptimizer<SteepestDescent, SteepestDescent> geo(*mockCalculator);
  geo.cellOptimizer.factor = 0.01;
  geo.geoOptimizer.optimizer.factor = 0.01;
  geo.check.maxIter = 1000;
  geo.geoMaxIterations = 900;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-8;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-8;
  geo.check.deltaValue = 1.0e-9;
  ASSERT_TRUE(GeometryOptimization::settingsMakeSense(geo));
  const ElementTypeCollection elements{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -0.05 * M_PI;
  positions(2, 1) = 0.95 * M_PI;
  AtomCollection atoms(elements, positions);
  auto pbc = PeriodicBoundaries(mockCalculator->settings().getString(SettingsNames::periodicBoundaries));
  // Def observer
  int counter = 0;
  auto func = [&](const int& /* param */, const double& /* value */, const Eigen::VectorXd& /* grad */) { counter++; };
  // Run
  geo.addObserver(func);
  const int nIter = geo.optimize(atoms, pbc, logger);
  // Check
  EXPECT_EQ(nIter, counter);
}

template<typename Optimizer>
void fixedAtomThrow() {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<UnitCellGeoOptMockCalculator>();
  UnitCellGeometryOptimizer<Optimizer, Optimizer> geo(*mockCalculator);
  geo.coordinateSystem = CoordinateSystem::Cartesian;
  geo.fixedAtoms = std::vector<int>{0, 1}; // Fix two atoms
  geo.check.maxIter = 50;

  const ElementTypeCollection elements{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -0.05 * M_PI;
  positions(2, 1) = 0.95 * M_PI;
  AtomCollection atoms(elements, positions);
  auto pbc = PeriodicBoundaries(mockCalculator->settings().getString(SettingsNames::periodicBoundaries));
  geo.optimize(atoms, pbc, logger);
}

TEST(UnitCellGeometryOptimizerTests, FixedAtomsThrow) {
  EXPECT_THROW(fixedAtomThrow<Bfgs>(), std::logic_error);
  EXPECT_THROW(fixedAtomThrow<Lbfgs>(), std::logic_error);
  EXPECT_THROW(fixedAtomThrow<SteepestDescent>(), std::logic_error);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
