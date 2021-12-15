/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/GeometryOptimization/GeometryOptimization.h"
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
#include <Utils/CalculatorBasics/TestCalculator.h>
#include <Utils/Constants.h>
#include <Utils/LennardJonesCalculator/LennardJonesCalculator.h>
#include <Utils/Settings.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

// Define mock calculator settings
class MockOptSettings : public Settings {
 public:
  MockOptSettings() : Settings("GeoOptMockCalculatorSettings") {
    UniversalSettings::DoubleDescriptor convergence_threshold("Energy convergence limit.");
    convergence_threshold.setDefaultValue(1e-12);
    this->_fields.push_back(SettingsNames::selfConsistenceCriterion, convergence_threshold);
    this->resetToDefaults();
  };
  ~MockOptSettings() override = default;
};

// Define a mock calculator
class GeoOptMockCalculator : public Core::Calculator {
 public:
  GeoOptMockCalculator() {
    settings_ = std::make_unique<MockOptSettings>();
  };
  ~GeoOptMockCalculator() override = default;
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
    return Utils::Property::Energy | Utils::Property::Gradients;
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

 private:
  AtomCollection structure_;
  Results r_;
  std::unique_ptr<Settings> settings_;
  Core::Calculator* cloneImpl() const final {
    return nullptr;
  }
};

/**
 * @class Scine::Utils::Tests::GeometryOptimizerTests
 * @brief Comprises tests for the class Scine::Utils::GeometryOptimizer.
 * @test
 */
TEST(GeometryOptimizerTests, SteepestDescent) {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<GeoOptMockCalculator>();
  GeometryOptimizer<SteepestDescent> geo(*mockCalculator);
  geo.coordinateSystem = CoordinateSystem::Cartesian;
  geo.check.maxIter = 50;
  geo.optimizer.factor = 0.01;
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
}

TEST(GeometryOptimizerTests, Bfgs) {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<GeoOptMockCalculator>();
  GeometryOptimizer<Bfgs> geo(*mockCalculator);
  geo.coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
  geo.check.maxIter = 20;
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
}

TEST(GeometryOptimizerTests, Bfgs2) {
  Core::Log logger = Core::Log::silent();
  TestCalculator testCalc;
  GeometryOptimizer<Bfgs> geo(testCalc);
  geo.coordinateSystem = CoordinateSystem::Cartesian;
  geo.check.maxIter = 200;
  geo.check.stepMaxCoeff = 1.0e-9;
  geo.check.stepRMS = 1.0e-10;
  geo.check.gradMaxCoeff = 1.0e-9;
  geo.check.gradRMS = 1.0e-10;
  geo.check.deltaValue = 1.0e-11;
  const ElementTypeCollection elements{ElementType::Br, ElementType::C, ElementType::H, ElementType::H, ElementType::H};
  constexpr double a2b = Constants::bohr_per_angstrom;
  PositionCollection positions = PositionCollection::Zero(5, 3);
  // clang-format off
  positions <<-3.6039746784*a2b,   +0.0000300140*a2b,   -0.1067055729*a2b,
              -1.6020084932*a2b,   -0.0000608235*a2b,   -0.0512927183*a2b,
              -1.2644784631*a2b,   +1.0265448353*a2b,   +0.1066552168*a2b,
              -1.2839109386*a2b,   -0.6416395584*a2b,   +0.7729317739*a2b,
              -1.2321296605*a2b,   -0.3848717385*a2b,   -1.0042935031*a2b;
  // clang-format on
  AtomCollection atoms(elements, positions);
  auto nIter = geo.optimize(atoms, logger);
  // Check results
  PositionCollection reference = PositionCollection::Zero(5, 3);
  // clang-format off
  reference << -7.1884413565e+00, +1.3427519106e-03, -2.1582324899e-01,
               -3.5853153517e+00, +6.7871512273e-05, -1.1227585886e-01,
               -2.0744585666e+00, +1.3260367874e+00, +1.2281489967e-01,
               -2.0999655137e+00, -8.2964494935e-01, +9.8383633867e-01,
               -2.0338472592e+00, -4.9779730464e-01, -1.3127867837e+00;
  // clang-format on
  positions = atoms.getPositions();
  for (unsigned int i = 0; i < reference.array().size(); i++) {
    EXPECT_NEAR(reference.data()[i], positions.data()[i], 1.0e-8);
  }
  ASSERT_EQ(nIter, 90);
}

TEST(GeometryOptimizerTests, Bfgs_NoDiis) {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<GeoOptMockCalculator>();
  GeometryOptimizer<Bfgs> geo(*mockCalculator);
  geo.coordinateSystem = CoordinateSystem::Cartesian;
  geo.optimizer.useGdiis = false;
  geo.check.maxIter = 30;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-8;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-8;
  geo.check.deltaValue = 1.0e-12;
  ASSERT_TRUE(GeometryOptimization::settingsMakeSense(geo));
  const ElementTypeCollection elements{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = 0.75 * M_PI;
  positions(2, 0) = 0.64 * M_PI;
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
}

TEST(GeometryOptimizerTests, Lbfgs) {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<GeoOptMockCalculator>();
  GeometryOptimizer<Lbfgs> geo(*mockCalculator);
  geo.coordinateSystem = CoordinateSystem::Cartesian;
  geo.optimizer.maxBacktracking = 20;
  geo.check.maxIter = 300;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-8;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-8;
  geo.check.deltaValue = 1.0e-12;
  ASSERT_TRUE(GeometryOptimization::settingsMakeSense(geo));
  const ElementTypeCollection elements{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = 0.75 * M_PI;
  positions(2, 0) = 0.64 * M_PI;
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
}

TEST(GeometryOptimizerTests, LbfgsNoLineSearch) {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<GeoOptMockCalculator>();
  GeometryOptimizer<Lbfgs> geo(*mockCalculator);
  geo.coordinateSystem = CoordinateSystem::Cartesian;
  geo.optimizer.linesearch = false;
  geo.check.maxIter = 20;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-9;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-9;
  geo.check.deltaValue = 1.0e-12;
  ASSERT_TRUE(GeometryOptimization::settingsMakeSense(geo));
  const ElementTypeCollection elements{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = 0.75 * M_PI;
  positions(2, 0) = 0.64 * M_PI;
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
}

TEST(GeometryOptimizerTests, NewtonRaphson) {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<GeoOptMockCalculator>();
  GeometryOptimizer<NewtonRaphson> geo(*mockCalculator);
  geo.check.maxIter = 50;
  geo.check.stepMaxCoeff = 1.0e-7;
  geo.check.stepRMS = 1.0e-8;
  geo.check.gradMaxCoeff = 1.0e-7;
  geo.check.gradRMS = 1.0e-8;
  geo.check.deltaValue = 1.0e-12;
  ASSERT_TRUE(GeometryOptimization::settingsMakeSense(geo));
  const ElementTypeCollection elements{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = -0.76 * M_PI;
  positions(2, 1) = 0.76 * M_PI;
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
}

TEST(GeometryOptimizerTests, TestObserver) {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<GeoOptMockCalculator>();
  GeometryOptimizer<SteepestDescent> geo(*mockCalculator);
  geo.optimizer.factor = 0.01;
  geo.check.maxIter = 63;
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
  // Def observer
  int counter = 0;
  auto func = [&](const int& /* param */, const double& /* value */, const Eigen::VectorXd& /* grad */) { counter++; };
  // Run
  geo.addObserver(func);
  const int nIter = geo.optimize(atoms, logger);
  // Check
  EXPECT_EQ(nIter, counter);
}

template<typename Optimizer>
bool fixedAtomConserved() {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<GeoOptMockCalculator>();
  GeometryOptimizer<Optimizer> geo(*mockCalculator);
  geo.coordinateSystem = CoordinateSystem::Cartesian;
  geo.fixedAtoms = std::vector<int>{0, 1}; // Fix two atoms
  geo.check.maxIter = 50;

  const ElementTypeCollection elements{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  const double zeroX = -0.05 * M_PI;
  const double twoY = 0.95 * M_PI;
  const double delta = 1e-2;
  positions(0, 0) = zeroX;
  positions(1, 2) = delta;
  positions(2, 1) = twoY - delta;
  AtomCollection atoms(elements, positions);

  try {
    geo.optimize(atoms, logger);
  }
  catch (std::exception& e) {
  }

  const auto& optimizedPositions = atoms.getPositions();

  return (optimizedPositions.row(0).isApprox(positions.row(0), 1e-8) &&
          optimizedPositions.row(1).isApprox(positions.row(1), 1e-8) &&
          !optimizedPositions.row(2).isApprox(positions.row(2), 1e-8));
}

TEST(GeometryOptimizerTests, FixedAtomsRespected) {
  // Optimize to minimum
  EXPECT_TRUE(fixedAtomConserved<SteepestDescent>());
  EXPECT_TRUE(fixedAtomConserved<NewtonRaphson>());
  EXPECT_TRUE(fixedAtomConserved<Lbfgs>());
  EXPECT_TRUE(fixedAtomConserved<Bfgs>());

  // Optimize to saddle point
  EXPECT_TRUE(fixedAtomConserved<Bofill>());
  EXPECT_TRUE(fixedAtomConserved<Dimer>());
  EXPECT_TRUE(fixedAtomConserved<EigenvectorFollowing>());
}

template<typename Optimizer>
bool fullyFixedTerminatesFast() {
  Core::Log logger = Core::Log::silent();
  std::shared_ptr<Core::Calculator> mockCalculator = std::make_shared<GeoOptMockCalculator>();
  GeometryOptimizer<Optimizer> geo(*mockCalculator);
  geo.coordinateSystem = CoordinateSystem::Cartesian;
  geo.fixedAtoms = std::vector<int>{0, 1, 2}; // Fix all atoms
  geo.check.maxIter = 100;

  const ElementTypeCollection elements{ElementType::H, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  const double zeroX = -0.05 * M_PI;
  const double twoY = 0.95 * M_PI;
  const double delta = 1e-2;
  positions(0, 0) = zeroX;
  positions(1, 2) = delta;
  positions(2, 1) = twoY - delta;
  AtomCollection atoms(elements, positions);

  int maxCycle = 0;
  geo.optimizer.addObserver([&maxCycle](const int& cycle, const double& /* value */,
                                        const Eigen::VectorXd& /* gradient */) { maxCycle = std::max(maxCycle, cycle); });

  try {
    geo.optimize(atoms, logger);
  }
  catch (std::exception& e) {
  }

  const auto& optimizedPositions = atoms.getPositions();

  const bool positionsUnchanged = (optimizedPositions.row(0).isApprox(positions.row(0), 1e-8) &&
                                   optimizedPositions.row(1).isApprox(positions.row(1), 1e-8) &&
                                   optimizedPositions.row(2).isApprox(positions.row(2), 1e-8));

  const bool fewCycles = maxCycle < 10;

  return positionsUnchanged && fewCycles;
}

TEST(GeometryOptimizerTests, FullyFixedTerminatesFast) {
  // Optimize to minimum
  EXPECT_TRUE(fullyFixedTerminatesFast<SteepestDescent>());
  EXPECT_TRUE(fullyFixedTerminatesFast<NewtonRaphson>());
  EXPECT_TRUE(fullyFixedTerminatesFast<Lbfgs>());
  EXPECT_TRUE(fullyFixedTerminatesFast<Bfgs>());
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
