/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/Interfaces/Calculator.h>
#include <Core/Log.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/MolecularDynamics/LeapFrogMD.h>
#include <Utils/MolecularDynamics/MolecularDynamics.h>
#include <Utils/MolecularDynamics/MolecularDynamicsSettings.h>
#include <Utils/Technical/CloneInterface.h>
#include <gmock/gmock.h>
#include <Eigen/Geometry>
#include <random>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

// Define a mock calculator
class MockCalculator : public CloneInterface<MockCalculator, Core::Calculator> {
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
    return PropertyList{};
  };
  PropertyList possibleProperties() const final {
    return Utils::Property::Energy | Utils::Property::Gradients;
  };
  const Results& calculate(std::string description) final {
    auto energy = 42.0;
    srand(42);
    Utils::GradientCollection gradients = Eigen::MatrixXd::Random(structure_.size(), 3);
    r_.set<Property::Energy>(energy);
    r_.set<Property::Gradients>(gradients);
    return r_;
  };
  std::string name() const final {
    return std::string("GeoOptMockCalculator");
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
    return nullptr;
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
 * @class AMolecularDynamicsTest MolecularDynamicsTest.cpp
 * @brief Comprises tests for Molecular Dynamics.
 * @test
 */
class AMolecularDynamicsTest : public Test {
 public:
  Core::Log log = Core::Log::silent();
};

TEST_F(AMolecularDynamicsTest, SettingsWorkCorrectly) {
  MockCalculator mockCalculator;
  MolecularDynamics md(mockCalculator);
  md.settings().modifyDouble(SettingsNames::timeStepInFemtoseconds, 3.0);
  md.settings().modifyInt(SettingsNames::numberOfMDSteps, 234);
  md.settings().modifyString(SettingsNames::integrationAlgorithm, "euler");
  md.settings().modifyBool(SettingsNames::temperatureBath, false);
  md.settings().modifyDouble(SettingsNames::relaxationTimeFactor, 15.0);
  md.settings().modifyBool(SettingsNames::saveVelocities, true);

  ASSERT_THAT(md.settings().getDouble(SettingsNames::timeStepInFemtoseconds), Eq(3.0));
  ASSERT_THAT(md.settings().getInt(SettingsNames::numberOfMDSteps), Eq(234));
  ASSERT_THAT(md.settings().getString(SettingsNames::integrationAlgorithm), Eq("euler"));
  ASSERT_FALSE(md.settings().getBool(SettingsNames::temperatureBath));
  ASSERT_THAT(md.settings().getDouble(SettingsNames::relaxationTimeFactor), Eq(15.0));
  ASSERT_TRUE(md.settings().getBool(SettingsNames::saveVelocities));
}

TEST_F(AMolecularDynamicsTest, mdSimulationIsPerformedCorrectlyWithLeapFrog) {
  MockCalculator mockCalculator;
  MolecularDynamics md(mockCalculator);

  md.settings().modifyInt(SettingsNames::numberOfMDSteps, 7);
  md.settings().modifyString(SettingsNames::integrationAlgorithm, OptionNames::leapFrogOption);
  md.settings().modifyBool(SettingsNames::saveVelocities, true);

  std::stringstream ss("3\n\n"
                       "O      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.4005871485    0.3100978455    0.0000000000\n"
                       "H     -0.4005871485    0.3100978455    0.0000000000\n");
  auto structure = Utils::XyzStreamHandler::read(ss);

  // Normal MD test
  md.performMDSimulation(structure, log);
  MolecularTrajectory mt = md.getMolecularTrajectory();

  ASSERT_THAT(mt.size(), Eq(8));
  for (const auto& energy : mt.getEnergies()) {
    ASSERT_THAT(energy, Eq(42.0));
  }
  for (const auto& pos : mt) {
    ASSERT_THAT(pos.rows(), Eq(structure.size()));
  }

  // Different temperature should lead to different trajectory
  md.settings().modifyDouble(SettingsNames::targetTemperature, 700);
  md.performMDSimulation(structure, log);
  MolecularTrajectory mtNew = md.getMolecularTrajectory();

  ASSERT_THAT(mtNew[0](1, 1), Eq(mt[0](1, 1)));
  for (int i = 1; i < mtNew.size(); ++i) {
    ASSERT_THAT(mtNew[i](1, 1), Ne(mt[i](1, 1)));
  }

  // Initial velocities should lead to different trajectory
  md.settings().modifyDouble(SettingsNames::targetTemperature, 300);
  md.setInitialVelocities(Eigen::MatrixXd::Random(structure.size(), 3));
  md.performMDSimulation(structure, log);
  MolecularTrajectory mtRandomStart = md.getMolecularTrajectory();

  ASSERT_THAT(mtRandomStart[0](1, 1), Eq(mt[0](1, 1)));
  for (int i = 1; i < mtRandomStart.size(); ++i) {
    ASSERT_THAT(mtRandomStart[i](1, 1), Ne(mt[i](1, 1)));
  }
  // Assert that velocities were saved as specified in the settings.
  ASSERT_THAT(md.getVelocities().size(), Eq(mtNew.size()));
}

TEST_F(AMolecularDynamicsTest, mdSimulationIsPerformedCorrectlyWithEuler) {
  MockCalculator mockCalculator;
  MolecularDynamics md(mockCalculator);

  md.settings().modifyInt(SettingsNames::numberOfMDSteps, 7);
  md.settings().modifyString(SettingsNames::integrationAlgorithm, OptionNames::eulerOption);

  std::stringstream ss("5\n\n"
                       "C      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.6287000000    0.6287000000    0.6287000000\n"
                       "H     -0.6287000000   -0.6287000000    0.6287000000\n"
                       "H     -0.6287000000    0.6287000000   -0.6287000000\n"
                       "H      0.6287000000   -0.6287000000   -0.6287000000\n");

  auto structure = Utils::XyzStreamHandler::read(ss);

  // Normal MD test
  md.performMDSimulation(structure, log);
  MolecularTrajectory mt = md.getMolecularTrajectory();

  ASSERT_THAT(mt.size(), Eq(8));
  for (const auto& energy : mt.getEnergies()) {
    ASSERT_THAT(energy, Eq(42.0));
  }
  for (const auto& pos : mt) {
    ASSERT_THAT(pos.rows(), Eq(structure.size()));
  }

  // Different temperature should lead to different trajectory
  md.settings().modifyDouble(SettingsNames::targetTemperature, 700);
  md.performMDSimulation(structure, log);
  MolecularTrajectory mtNew = md.getMolecularTrajectory();

  ASSERT_THAT(mtNew[0](1, 1), Eq(mt[0](1, 1)));
  // Start at 2 because 1st structure is also equal for Euler algorithm.
  for (int i = 2; i < mtNew.size(); ++i) {
    ASSERT_THAT(mtNew[i](1, 1), Ne(mt[i](1, 1)));
  }

  // Initial velocities should lead to different trajectory
  md.settings().modifyDouble(SettingsNames::targetTemperature, 300);
  md.setInitialVelocities(Eigen::MatrixXd::Random(structure.size(), 3));
  md.performMDSimulation(structure, log);
  MolecularTrajectory mtRandomStart = md.getMolecularTrajectory();

  ASSERT_THAT(mtRandomStart[0](1, 1), Eq(mt[0](1, 1)));
  for (int i = 1; i < mtRandomStart.size(); ++i) {
    ASSERT_THAT(mtRandomStart[i](1, 1), Ne(mt[i](1, 1)));
  }
  // Should be empty since the default setting is to not record the velocities.
  ASSERT_TRUE(md.getVelocities().empty());
}

TEST_F(AMolecularDynamicsTest, mdSimulationIsPerformedCorrectlyWithVelocityVerlet) {
  MockCalculator mockCalculator;
  MolecularDynamics md(mockCalculator);

  md.settings().modifyInt(SettingsNames::numberOfMDSteps, 7);
  md.settings().modifyString(SettingsNames::integrationAlgorithm, OptionNames::velocityVerletOption);

  std::stringstream ss("2\n\n"
                       "H     0.0000000000    0.0000000000   -0.0000000000\n"
                       "H     2.0000000000    0.0000000000    0.0000000000\n");

  auto structure = Utils::XyzStreamHandler::read(ss);

  // Normal MD test
  md.performMDSimulation(structure, log);
  MolecularTrajectory mt = md.getMolecularTrajectory();

  ASSERT_THAT(mt.size(), Eq(8));
  for (const auto& energy : mt.getEnergies()) {
    ASSERT_THAT(energy, Eq(42.0));
  }
  for (const auto& pos : mt) {
    ASSERT_THAT(pos.rows(), Eq(structure.size()));
  }

  // Different temperature should lead to different trajectory
  md.settings().modifyDouble(SettingsNames::targetTemperature, 700);
  md.performMDSimulation(structure, log);
  MolecularTrajectory mtNew = md.getMolecularTrajectory();

  ASSERT_THAT(mtNew[0](1, 1), Eq(mt[0](1, 1)));
  // Start at 2 because 1st structure is also equal for Velocity Verlet algorithm.
  for (int i = 2; i < mtNew.size(); ++i) {
    ASSERT_THAT(mtNew[i](1, 1), Ne(mt[i](1, 1)));
  }

  // Initial velocities should lead to different trajectory
  md.settings().modifyDouble(SettingsNames::targetTemperature, 300);
  md.setInitialVelocities(Eigen::MatrixXd::Random(structure.size(), 3));
  md.performMDSimulation(structure, log);
  MolecularTrajectory mtRandomStart = md.getMolecularTrajectory();

  ASSERT_THAT(mtRandomStart[0](1, 1), Eq(mt[0](1, 1)));
  for (int i = 1; i < mtRandomStart.size(); ++i) {
    ASSERT_THAT(mtRandomStart[i](1, 1), Ne(mt[i](1, 1)));
  }
  // Should be empty since the default setting is to not record the velocities.
  ASSERT_TRUE(md.getVelocities().empty());
}

TEST_F(AMolecularDynamicsTest, InitialVelocitiesGeneratedAndCOMComponentRemoved) {
  LeapFrogMD mdEngine;
  mdEngine.setSeed(42);
  ElementTypeCollection elements = {ElementType::H};
  mdEngine.setElementTypes(elements);

  PositionCollection h2Positions;
  h2Positions.resize(2, 3);
  h2Positions << 1.6, 0, 0, 0, 0, 0;
  mdEngine.removeCenterOfMassLinearMomentum(h2Positions.row(0));
  EXPECT_THAT(mdEngine.getVelocities().row(0)(0), DoubleNear(0, 1e-9));
  EXPECT_THAT(mdEngine.getVelocities().row(0)(1), DoubleNear(0, 1e-9));
  EXPECT_THAT(mdEngine.getVelocities().row(0)(2), DoubleNear(0, 1e-9));

  ElementTypeCollection elements2 = {ElementType::C, ElementType::H};
  mdEngine.setElementTypes(elements2);
  mdEngine.removeCenterOfMassLinearMomentum(h2Positions);
  mdEngine.removeCenterOfMassAngularMomentum(h2Positions);
  Eigen::MatrixX3d velocities2 = mdEngine.getVelocities();
  auto masses = Utils::Geometry::getMasses(elements2);
  auto centerOfMass = Utils::Geometry::getCenterOfMass(h2Positions, masses);
  Eigen::MatrixX3d linearMomenta =
      velocities2.array().colwise() * Eigen::Map<const Eigen::ArrayXd>(masses.data(), masses.size());
  Eigen::RowVector3d linearMomentum = linearMomenta.colwise().sum();
  EXPECT_THAT(linearMomentum.x(), DoubleNear(0, 1e-9));
  EXPECT_THAT(linearMomentum.y(), DoubleNear(0, 1e-9));
  EXPECT_THAT(linearMomentum.z(), DoubleNear(0, 1e-9));

  Eigen::RowVector3d angularMomentum = Eigen::RowVector3d::Zero();
  for (int i = 0; i < 2; ++i) {
    angularMomentum += (h2Positions.row(i) - centerOfMass).cross(linearMomenta.row(i));
  }

  EXPECT_THAT(angularMomentum.x(), DoubleNear(0, 1e-9));
  EXPECT_THAT(angularMomentum.y(), DoubleNear(0, 1e-9));
  EXPECT_THAT(angularMomentum.z(), DoubleNear(0, 1e-9));

  PositionCollection arbitraryPositions = 10. * Eigen::MatrixX3d::Random(5, 3);
  ElementTypeCollection arbitraryElements = {ElementType::H, ElementType::C, ElementType::P, ElementType::Al, ElementType::Ar};

  mdEngine.setElementTypes(arbitraryElements);
  mdEngine.removeCenterOfMassLinearMomentum(arbitraryPositions);
  mdEngine.removeCenterOfMassAngularMomentum(arbitraryPositions);

  Eigen::MatrixX3d v3 = mdEngine.getVelocities();
  auto masses3 = Utils::Geometry::getMasses(arbitraryElements);
  auto centerOfMass3 = Utils::Geometry::getCenterOfMass(arbitraryPositions, masses3);
  Eigen::MatrixX3d linearMomenta3 = v3.array().colwise() * Eigen::Map<const Eigen::ArrayXd>(masses3.data(), masses3.size());
  Eigen::RowVector3d linearMomentum3 = linearMomenta3.colwise().sum();
  EXPECT_THAT(linearMomentum3.x(), DoubleNear(0, 1e-9));
  EXPECT_THAT(linearMomentum3.y(), DoubleNear(0, 1e-9));
  EXPECT_THAT(linearMomentum3.z(), DoubleNear(0, 1e-9));
  Eigen::RowVector3d angularMomentum3 = Eigen::RowVector3d::Zero();

  for (int i = 0; i < 5; ++i) {
    angularMomentum3 += (arbitraryPositions.row(i) - centerOfMass3).cross(linearMomenta3.row(i));
  }

  EXPECT_THAT(angularMomentum3.x(), DoubleNear(0, 1e-9));
  EXPECT_THAT(angularMomentum3.y(), DoubleNear(0, 1e-9));
  EXPECT_THAT(angularMomentum3.z(), DoubleNear(0, 1e-9));
}

TEST_F(AMolecularDynamicsTest, InitialVelocitiesAreGeneratedAtCorrectTemperature) {
  LeapFrogMD mdEngine;
  mdEngine.setTargetTemperatureInKelvin(298.15);
  ElementTypeCollection arbitraryElements = {
      ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,
      ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,
      ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,
      ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al,
      ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar,
      ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,
      ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,
      ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,
      ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al,
      ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar,
      ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,
      ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,
      ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,
      ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al,
      ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar,
      ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,
      ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,
      ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,
      ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al,
      ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar,
      ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,
      ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,
      ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,
      ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al,
      ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar,
      ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,
      ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,
      ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,
      ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al,
      ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar,
      ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,
      ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,
      ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,
      ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al,
      ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar,
      ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,
      ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,
      ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,
      ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al,
      ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar,
      ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,
      ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,
      ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,
      ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al,
      ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar,
      ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,
      ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,
      ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,
      ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al,
      ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar,
      ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,
      ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,
      ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,
      ElementType::Al, ElementType::Ar};
  PositionCollection arbitraryPositions = 10. * Eigen::MatrixX3d::Random(arbitraryElements.size(), 3);

  std::vector<int> seeds(1000);
  std::vector<double> temperatures;
  temperatures.reserve(1000);
  std::iota(seeds.begin(), seeds.end(), 123958);

  auto masses = Utils::Geometry::getMasses(arbitraryElements);

  for (int seed : seeds) {
    mdEngine.setSeed(seed);
    mdEngine.setElementTypes(arbitraryElements);

    auto velocities = mdEngine.getVelocities();
    double temperature =
        (velocities.rowwise().squaredNorm().array() * Eigen::Map<const Eigen::ArrayXd>(masses.data(), masses.size())).sum() /
        (3. * masses.size());
    temperatures.push_back(temperature);
  }

  double meanTemperature = std::accumulate(temperatures.begin(), temperatures.end(), 0.) / temperatures.size();
  ASSERT_THAT(meanTemperature * 3.1577464e5, DoubleNear(298.15, 1));
}
} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
