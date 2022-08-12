/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Tested File */
#include "Utils/GeometryOptimization/NtOptimizer2.h"
/* Internal */
#include "Core/Interfaces/Calculator.h"
#include "Utils/Bonds/BondDetector.h"
#include "Utils/CalculatorBasics/PropertyList.h"
#include "Utils/CalculatorBasics/Results.h"
#include "Utils/CalculatorBasics/TestCalculator.h"
#include "Utils/IO/MolecularTrajectoryIO.h"
#include "Utils/Math/QuaternionFit.h"
#include "Utils/Settings.h"
#include "Utils/Technical/CloneInterface.h"
/* External */
#include <gmock/gmock.h>
#include <boost/dll/runtime_symbol_info.hpp>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

// Define a mock calculator
class NtOpt2MockCalculator : public CloneInterface<NtOpt2MockCalculator, Core::Calculator> {
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
    return Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::BondOrderMatrix;
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
    // Bond Orders
    auto bos = BondDetector::detectBonds(structure_);

    r_ = Results();
    r_.set<Property::SuccessfulCalculation>(true);
    r_.set<Property::Energy>(e);
    r_.set<Property::Gradients>(g);
    r_.set<Property::BondOrderMatrix>(bos);
    return r_;
  };
  std::string name() const final {
    return std::string("NtOpt2MockCalculator");
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
 * @class Scine::Utils::Tests::NtOptimizer2Tests
 * @brief Comprises tests for the class Scine::Utils::NtOptimizer2.
 * @test
 */
TEST(NtOptimizer2Tests, DefaultEmptyLists) {
  Core::Log logger = Core::Log::silent();
  NtOpt2MockCalculator mockCalculator;
  NtOptimizer2 nt(mockCalculator);
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(2, 3);
  positions(0, 1) = 11.0;
  AtomCollection atoms(elements, positions);
  ASSERT_THROW(nt.optimize(atoms, logger), std::logic_error);
}

inline void runSN2Test(NtOptimizer2& nt, Core::Log& logger) {
  nt.coordinateSystem = CoordinateSystem::Cartesian;
  nt.associationList = {0, 1};
  nt.dissociationList = {1, 5};
  nt.useMicroCycles = true;
  nt.totalForceNorm = 0.1;
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
  std::vector<int> referenceReactiveAtomsList = {0, 1, 5};
  std::vector<int> referenceAxesList = {0, 1, 1, 5};

  for (unsigned int i = 0; i < referenceReactiveAtomsList.size(); i++) {
    ASSERT_EQ(nt.getReactiveAtomsList()[i], referenceReactiveAtomsList[i]);
  }
  for (unsigned int i = 0; i < referenceAxesList.size(); i++) {
    ASSERT_EQ(nt.getConstraintsMap()[1][i], referenceAxesList[i]);
  };
  ASSERT_EQ(nIter, 37);
  // Check results
  PositionCollection reference = PositionCollection::Zero(6, 3);
  Eigen::VectorXd referenceBonds = Eigen::VectorXd::Zero(4, 1);
  // clang-format off
  referenceBonds << 4.5981931756e+00, //Cl-C
                    9.8883569663e+00, //Cl-Br
                    2.0063532273e+00, //C-H
                    5.2901692863e+00; //C-Br
  // clang-format on
  PositionCollection positions = atoms.getPositions();

  EXPECT_NEAR(referenceBonds.row(0).value(), (positions.row(0) - positions.row(1)).norm(), 1e-3); // Cl-C
  EXPECT_NEAR(referenceBonds.row(1).value(), (positions.row(0) - positions.row(5)).norm(), 1e-3); // Cl-Br
  EXPECT_NEAR(referenceBonds.row(2).value(), (positions.row(1) - positions.row(4)).norm(), 1e-3); // C-H
  EXPECT_NEAR(referenceBonds.row(3).value(), (positions.row(1) - positions.row(5)).norm(), 1e-3); // C-Br
}

TEST(NtOptimizer2Tests, SN2) {
  Core::Log logger = Core::Log::silent();
  TestCalculator testCalc;
  NtOptimizer2 nt(testCalc);
  runSN2Test(nt, logger);
}

TEST(NtOptimizer2Tests, SN2FirstMaximum) {
  Core::Log logger = Core::Log::silent();
  TestCalculator testCalc;
  NtOptimizer2 nt(testCalc);
  nt.extractionCriterion = NtOptimizer2::ntExtractFirst;
  runSN2Test(nt, logger);
}

TEST(NtOptimizer2Tests, SN2HighestMaximum) {
  Core::Log logger = Core::Log::silent();
  TestCalculator testCalc;
  NtOptimizer2 nt(testCalc);
  nt.extractionCriterion = NtOptimizer2::ntExtractHighest;
  runSN2Test(nt, logger);
}

TEST(NtOptimizer2Tests, SingleStepFailure) {
  Core::Log logger = Core::Log::silent();
  NtOpt2MockCalculator mockCalculator;
  NtOptimizer2 nt(mockCalculator);
  nt.associationList = {0, 1};
  nt.totalForceNorm = 1.0;
  nt.check.maxIter = 10;
  nt.optimizer.factor = 0.02;
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(2, 3);
  // Enforce a scan that stops after the initial evaluation.
  positions(0, 1) = 0.1;
  AtomCollection atoms(elements, positions);
  // Expect failure due to missing TS (trajectory has only one point)
  ASSERT_THROW(nt.optimize(atoms, logger), std::runtime_error);
}

TEST(NtOptimizer2Tests, ConnectedNucleiWorks) {
  Core::Log logger = Core::Log::silent();
  NtOpt2MockCalculator mockCalculator;
  mockCalculator.setLog(logger);
  auto pathToResources = boost::dll::program_location().parent_path();
  pathToResources /= "Resources";
  const auto trajectoryFile = (pathToResources / "dissociation.trj.xyz").string();
  auto trajectory = Utils::MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::xyz, trajectoryFile);
  // testing connected state
  auto atomsFirst = Utils::AtomCollection(trajectory.getElementTypes(), trajectory.front());
  mockCalculator.setStructure(atomsFirst);
  auto result = mockCalculator.calculate();
  ASSERT_TRUE(result.has<Utils::Property::BondOrderMatrix>());
  const auto boFirst = result.get<Utils::Property::BondOrderMatrix>();
  std::vector<int> indices = {0, 2, 4};
  auto molecules = NtUtils::connectedNuclei(indices, boFirst);
  ASSERT_THAT(molecules.size(), 1);
  ASSERT_THAT(molecules[0], indices);
  // testing unconnected state
  auto atomsLast = Utils::AtomCollection(trajectory.getElementTypes(), trajectory.back());
  mockCalculator.setStructure(atomsLast);
  result = mockCalculator.calculate();
  ASSERT_TRUE(result.has<Utils::Property::BondOrderMatrix>());
  const auto boLast = result.get<Utils::Property::BondOrderMatrix>();
  indices = {1, 9};
  molecules = NtUtils::connectedNuclei(indices, boLast);
  ASSERT_THAT(molecules.size(), 2);
  ASSERT_TRUE((molecules[0] == std::vector<int>{indices[0]} && molecules[1] == std::vector<int>{indices[1]}) ||
              (molecules[1] == std::vector<int>{indices[0]} && molecules[0] == std::vector<int>{indices[1]}));
}

TEST(NtOptimizer2Tests, Center2CenterWorks) {
  PositionCollection pos = PositionCollection::Zero(4, 3);
  // clang-format off
  pos <<  -1.0,  0.0,  0.0,
           1.0,  2.0, 10.0,
           1.0, -2.0, 10.0,
           4.2,  3.1, -0.9;
  // clang-format on
  Displacement expected;
  expected << -2.0, 0.0, -10.0;
  std::vector<int> lhs = {0};
  std::vector<int> rhs = {1, 2};
  auto c2c = NtUtils::centerToCenterVector(pos, lhs, rhs);
  ASSERT_TRUE(expected.isApprox(c2c));
}

TEST(NtOptimizer2Tests, SmallestCovRadiusWorks) {
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
  std::vector<int> i = {0, 1, 5};
  std::vector<int> j = {1, 2, 3};
  ASSERT_DOUBLE_EQ(NtUtils::smallestCovalentRadius(atoms, i), ElementInfo::covalentRadius(ElementType::C));
  ASSERT_DOUBLE_EQ(NtUtils::smallestCovalentRadius(atoms, j), ElementInfo::covalentRadius(ElementType::H));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
