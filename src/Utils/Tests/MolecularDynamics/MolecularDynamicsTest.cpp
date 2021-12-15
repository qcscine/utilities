/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/CalculatorWithReference.h>
#include <Core/Log.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/CalculatorBasics/TestCalculator.h>
#include <Utils/Constants.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/MolecularDynamics/LeapFrogMD.h>
#include <Utils/MolecularDynamics/MolecularDynamics.h>
#include <Utils/MolecularDynamics/MolecularDynamicsSettings.h>
#include <Utils/Technical/CloneInterface.h>
#include <gmock/gmock.h>
#include <Eigen/Geometry>
#include <functional>
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
    return PropertyList{};
  };
  PropertyList possibleProperties() const final {
    return Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::BondOrderMatrix |
           Utils::Property::AtomicCharges;
  };

  const Results& calculate(std::string /*description*/) final {
    auto energy = 42.0;
    srand(42);
    Utils::GradientCollection gradients = Eigen::MatrixXd::Random(structure_.size(), 3);
    Utils::BondOrderCollection bondOrders(structure_.size());
    std::vector<double> atomicCharges(structure_.size(), 42.0);
    r_.set<Property::Energy>(energy);
    r_.set<Property::Gradients>(gradients);
    r_.set<Property::BondOrderMatrix>(bondOrders);
    r_.set<Property::AtomicCharges>(atomicCharges);
    r_.set<Property::SuccessfulCalculation>(true);
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
  void loadState(std::shared_ptr<Core::State> /* state */) final {
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
  Settings settings_ = Settings("dummy");
};

// Define a mock calculator with reference calculator
class MockCalculatorWithReference final : public Core::CalculatorWithReference {
 public:
  const Results& calculate() final {
    referenceCalculation();
    if (!referenceCalculator_) {
      throw std::runtime_error("Mock calculator with reference needs reference calculator");
    }

    auto results = referenceCalculator_->results();
    srand(42);
    double energy = results.get<Property::Energy>() * 2.0;
    Utils::GradientCollection gradients = results.get<Property::Gradients>();
    Utils::BondOrderCollection bondOrders = results.get<Property::BondOrderMatrix>();
    bondOrders.setOrder(0, 1, 3.0);
    std::vector<double> atomicCharges = results.get<Property::AtomicCharges>();
    std::transform(atomicCharges.begin(), atomicCharges.end(), atomicCharges.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, 0.5));

    r_.set<Property::Energy>(energy);
    r_.set<Property::Gradients>(gradients);
    r_.set<Property::BondOrderMatrix>(bondOrders);
    r_.set<Property::AtomicCharges>(atomicCharges);
    r_.set<Property::SuccessfulCalculation>(true);
    return r_;
  };

  void setReferenceCalculator(std::shared_ptr<Core::Calculator> referenceCalculator) {
    referenceCalculator_ = std::dynamic_pointer_cast<Core::Calculator>(referenceCalculator);
  }

  void referenceCalculation() {
    referenceCalculator_->calculate();
  }

  Core::Calculator& getReferenceCalculator() {
    return *referenceCalculator_;
  }

  const Core::Calculator& getReferenceCalculator() const {
    return *referenceCalculator_;
  }

  std::string name() const final {
    return std::string("MockCalculatorWithReference");
  };
  const Settings& settings() const final {
    return *settings_;
  }
  Settings& settings() final {
    return *settings_;
  }
  void applySettings() {
  }

  Utils::Results& results() final {
    return r_;
  }
  const Utils::Results& results() const final {
    return r_;
  }

 private:
  Results r_;
  std::shared_ptr<Core::Calculator> referenceCalculator_;
  std::unique_ptr<Settings> settings_ = std::make_unique<Settings>(Settings("dummy"));
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
  md.settings().modifyDouble(SettingsNames::generationTemperature, 200.0);
  md.settings().modifyInt(SettingsNames::generationSeed, 3);
  md.settings().modifyString(SettingsNames::thermostatAlgorithm, "berendsen");
  md.settings().modifyDouble(SettingsNames::targetTemperature, 500.0);
  md.settings().modifyDouble(SettingsNames::temperatureCouplingTime, 15.0);
  md.settings().modifyInt(SettingsNames::stochasticDynamicsSeed, 4);
  md.settings().modifyBool(SettingsNames::saveVelocities, true);
  md.settings().modifyBool(SettingsNames::saveTemperatures, true);
  md.settings().modifyBool(SettingsNames::requireCharges, true);
  md.settings().modifyBool(SettingsNames::requireBondOrders, true);

  ASSERT_THAT(md.settings().getDouble(SettingsNames::timeStepInFemtoseconds), Eq(3.0));
  ASSERT_THAT(md.settings().getInt(SettingsNames::numberOfMDSteps), Eq(234));
  ASSERT_THAT(md.settings().getString(SettingsNames::integrationAlgorithm), Eq("euler"));
  ASSERT_THAT(md.settings().getDouble(SettingsNames::generationTemperature), Eq(200.0));
  ASSERT_THAT(md.settings().getInt(SettingsNames::generationSeed), Eq(3));
  ASSERT_THAT(md.settings().getString(SettingsNames::thermostatAlgorithm), Eq("berendsen"));
  ASSERT_THAT(md.settings().getDouble(SettingsNames::targetTemperature), Eq(500.0));
  ASSERT_THAT(md.settings().getDouble(SettingsNames::temperatureCouplingTime), Eq(15.0));
  ASSERT_THAT(md.settings().getInt(SettingsNames::stochasticDynamicsSeed), Eq(4));
  ASSERT_TRUE(md.settings().getBool(SettingsNames::saveVelocities));
  ASSERT_TRUE(md.settings().getBool(SettingsNames::saveTemperatures));
  ASSERT_TRUE(md.settings().getBool(SettingsNames::requireCharges));
  ASSERT_TRUE(md.settings().getBool(SettingsNames::requireBondOrders));
}

TEST_F(AMolecularDynamicsTest, mdSimulationIsPerformedCorrectlyWithLeapFrog) {
  MockCalculator mockCalculator;
  MolecularDynamics md(mockCalculator);

  md.settings().modifyInt(SettingsNames::numberOfMDSteps, 7);
  md.settings().modifyString(SettingsNames::integrationAlgorithm, OptionNames::leapFrogOption);
  md.settings().modifyString(SettingsNames::thermostatAlgorithm, "berendsen");
  md.settings().modifyBool(SettingsNames::saveVelocities, true);
  md.settings().modifyBool(SettingsNames::saveTemperatures, true);

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
  int counter = 0;
  for (const auto& pos : mt) {
    ASSERT_THAT(pos.rows(), Eq(structure.size()));
    if (counter == 0)
      ASSERT_TRUE(pos.isApprox(structure.getPositions()));
    else
      ASSERT_FALSE(pos.isApprox(structure.getPositions()));
    counter++;
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
  auto velocities = md.getVelocities();
  ASSERT_THAT(velocities.size(), Eq(mtNew.size()));
  // Check that final velocities can be retrieved correctly
  EXPECT_TRUE(md.getFinalVelocities().isApprox(velocities.back()));
  // Assert that temperatures were saved as specified in the settings.
  auto temperatures = md.getTemperatures();
  ASSERT_THAT(temperatures.size(), Eq(mtNew.size()));
  // Assert that stored temperatures fit with stored velocities
  auto masses = Utils::Geometry::Properties::getMasses(structure.getElements());
  double velocityTemperature;
  for (unsigned int i = 0; i < velocities.size(); i++) {
    velocityTemperature =
        (velocities[i].rowwise().squaredNorm().array() * Eigen::Map<const Eigen::ArrayXd>(masses.data(), masses.size())).sum() /
        (3. * masses.size());
    // Convert to Kelvin
    velocityTemperature *= (Utils::Constants::joule_per_hartree / Utils::Constants::boltzmannConstant);
    EXPECT_THAT(temperatures[i], DoubleNear(velocityTemperature, 1e-9));
  }
}

TEST_F(AMolecularDynamicsTest, mdSimulationIsPerformedCorrectlyWithEuler) {
  MockCalculator mockCalculator;
  MolecularDynamics md(mockCalculator);

  md.settings().modifyInt(SettingsNames::numberOfMDSteps, 7);
  md.settings().modifyString(SettingsNames::integrationAlgorithm, OptionNames::eulerOption);
  md.settings().modifyString(SettingsNames::thermostatAlgorithm, "berendsen");

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
  int counter = 0;
  for (const auto& pos : mt) {
    ASSERT_THAT(pos.rows(), Eq(structure.size()));
    if (counter == 0)
      ASSERT_TRUE(pos.isApprox(structure.getPositions()));
    else
      ASSERT_FALSE(pos.isApprox(structure.getPositions()));
    counter++;
  }

  // Different temperature should lead to different trajectory
  md.settings().modifyDouble(SettingsNames::targetTemperature, 700);
  md.performMDSimulation(structure, log);
  MolecularTrajectory mtNew = md.getMolecularTrajectory();

  // Final velocities should always be obtainable
  EXPECT_THAT(md.getFinalVelocities().rows(), Eq(structure.size()));

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
  // Should be empty since the default setting is to not record the temperatures.
  ASSERT_TRUE(md.getTemperatures().empty());
}

TEST_F(AMolecularDynamicsTest, mdSimulationIsPerformedCorrectlyWithVelocityVerlet) {
  MockCalculator mockCalculator;
  MolecularDynamics md(mockCalculator);

  md.settings().modifyInt(SettingsNames::numberOfMDSteps, 7);
  md.settings().modifyString(SettingsNames::integrationAlgorithm, OptionNames::velocityVerletOption);
  md.settings().modifyString(SettingsNames::thermostatAlgorithm, "berendsen");

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
  // Should be empty since the default setting is to not record the temperatures.
  ASSERT_TRUE(md.getTemperatures().empty());
}

TEST_F(AMolecularDynamicsTest, mdSimulationIsPerformedCorrectlyWithStochasticDynamics) {
  TestCalculator testCalculator;
  MolecularDynamics md(testCalculator);

  md.settings().modifyInt(SettingsNames::numberOfMDSteps, 200);
  md.settings().modifyDouble(SettingsNames::timeStepInFemtoseconds, 0.5);
  // Use small tau to observe rapid equilibration
  md.settings().modifyDouble(SettingsNames::temperatureCouplingTime, 0.5);
  md.settings().modifyString(SettingsNames::integrationAlgorithm, OptionNames::stochasticDynamicsOption);
  md.settings().modifyDouble(SettingsNames::generationTemperature, 700.0);
  md.settings().modifyDouble(SettingsNames::targetTemperature, 250.0);
  md.settings().modifyBool(SettingsNames::saveTemperatures, true);

  // 10 water molecules
  std::stringstream ss("30\n\n"
                       "O      -1.09205000      -0.67155000       1.79165000\n"
                       "H      -0.72504000      -1.14569000       2.53778000\n"
                       "H      -1.28155000       0.20261000       2.13248000\n"
                       "O       2.16941000      -1.40715000       1.63986000\n"
                       "H       1.43752000      -0.89847000       1.98887000\n"
                       "H       1.77185000      -2.22415000       1.33874000\n"
                       "O       2.57322000      -1.03407000      -1.72070000\n"
                       "H       2.49262000      -1.94861000      -1.99153000\n"
                       "H       2.15966000      -0.53736000      -2.42673000\n"
                       "O       0.14592000       3.41955000       0.16725000\n"
                       "H       0.07675000       2.47854000       0.00616000\n"
                       "H       0.90791000       3.50969000       0.73950000\n"
                       "O      -3.21061000      -0.00178000      -1.69577000\n"
                       "H      -3.61492000      -0.71405000      -2.19118000\n"
                       "H      -3.12430000      -0.34645000      -0.80696000\n"
                       "O      -2.85229000      -2.43221000       0.18840000\n"
                       "H      -2.34446000      -1.93412000       0.82890000\n"
                       "H      -2.32096000      -3.20869000       0.01234000\n"
                       "O       1.17644000      -3.77341000       0.39414000\n"
                       "H       1.55292000      -4.64622000       0.28151000\n"
                       "H       1.07839000      -3.43754000      -0.49682000\n"
                       "O       3.03950000       2.83166000      -1.01808000\n"
                       "H       3.39188000       2.23448000      -1.67794000\n"
                       "H       2.15573000       2.50609000      -0.84726000\n"
                       "O      -1.57145000       0.59603000      -3.99850000\n"
                       "H      -1.11518000       1.40656000      -3.77249000\n"
                       "H      -1.92263000       0.27835000      -3.16666000\n"
                       "O      -0.21681000      -2.15490000       3.90836000\n"
                       "H       0.28172000      -1.68032000       4.57354000\n"
                       "H      -0.07451000      -3.07855000       4.11536000\n");

  auto structure = Utils::XyzStreamHandler::read(ss);

  // Normal MD test
  md.performMDSimulation(structure, log);
  MolecularTrajectory mt = md.getMolecularTrajectory();

  // Formal tests
  ASSERT_THAT(mt.size(), Eq(201));
  for (const auto& pos : mt) {
    EXPECT_THAT(pos.rows(), Eq(structure.size()));
  }
  std::vector<double> temperatures = md.getTemperatures();
  ASSERT_THAT(temperatures.size(), Eq(201));
  EXPECT_TRUE(md.getVelocities().empty());
  // Temperatures should have been reduced to ca. 250
  double average = std::accumulate(temperatures.begin() + 20, temperatures.end(), 0.0) / 181.0;
  // Due to small system size and small number of steps notable deviations are to be expected
  EXPECT_THAT(average, DoubleNear(250, 50));
}

TEST_F(AMolecularDynamicsTest, InitialVelocitiesGeneratedAndCOMComponentRemoved) {
  LeapFrogMD mdEngine;
  mdEngine.setSeed(42);
  ElementTypeCollection elements = {ElementType::H};
  mdEngine.setElementTypes(elements);

  PositionCollection h2Positions;
  h2Positions.resize(2, 3);
  h2Positions << 1.6, 0, 0, 0, 0, 0;
  mdEngine.removeCenterOfMassLinearMomentum();
  EXPECT_THAT(mdEngine.getVelocities().row(0)(0), DoubleNear(0, 1e-9));
  EXPECT_THAT(mdEngine.getVelocities().row(0)(1), DoubleNear(0, 1e-9));
  EXPECT_THAT(mdEngine.getVelocities().row(0)(2), DoubleNear(0, 1e-9));

  ElementTypeCollection elements2 = {ElementType::C, ElementType::H};
  mdEngine.setElementTypes(elements2);
  mdEngine.removeCenterOfMassLinearMomentum();
  mdEngine.removeCenterOfMassAngularMomentum(h2Positions);
  Eigen::MatrixX3d velocities2 = mdEngine.getVelocities();
  auto masses = Utils::Geometry::Properties::getMasses(elements2);
  auto centerOfMass = Utils::Geometry::Properties::getCenterOfMass(h2Positions, masses);
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
  mdEngine.removeCenterOfMassLinearMomentum();
  mdEngine.removeCenterOfMassAngularMomentum(arbitraryPositions);

  Eigen::MatrixX3d v3 = mdEngine.getVelocities();
  auto masses3 = Utils::Geometry::Properties::getMasses(arbitraryElements);
  auto centerOfMass3 = Utils::Geometry::Properties::getCenterOfMass(arbitraryPositions, masses3);
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

  auto masses = Utils::Geometry::Properties::getMasses(arbitraryElements);

  for (int seed : seeds) {
    mdEngine.setSeed(seed);
    mdEngine.setElementTypes(arbitraryElements);
    mdEngine.setGenerationTemperatureInKelvin(400);
    mdEngine.sampleVelocitiesFromBoltzmannDistribution();

    auto velocities = mdEngine.getVelocities();
    double temperature =
        (velocities.rowwise().squaredNorm().array() * Eigen::Map<const Eigen::ArrayXd>(masses.data(), masses.size())).sum() /
        (3. * masses.size());
    temperatures.push_back(temperature);
  }

  double meanTemperature = std::accumulate(temperatures.begin(), temperatures.end(), 0.) / temperatures.size();
  ASSERT_THAT(meanTemperature * 3.1577464e5, DoubleNear(400., 1));
}

TEST_F(AMolecularDynamicsTest, mdSimulationRespectsInitialVelocities) {
  MockCalculator mockCalculator;
  MolecularDynamics md(mockCalculator);

  md.settings().modifyInt(SettingsNames::numberOfMDSteps, 1);
  md.settings().modifyBool(SettingsNames::saveVelocities, true);
  md.settings().modifyBool(SettingsNames::saveTemperatures, true);

  ElementTypeCollection arbitraryElements = {
      ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::N,
      ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,
      ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,
      ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al,
      ElementType::Ne, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar,
      ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Ar, ElementType::H,
      ElementType::C,  ElementType::P,  ElementType::Al, ElementType::Fe, ElementType::H,  ElementType::C,
      ElementType::P,  ElementType::Al, ElementType::F,  ElementType::H,  ElementType::C,  ElementType::P,
      ElementType::Al, ElementType::Ar, ElementType::H,  ElementType::C,  ElementType::P,  ElementType::Al};

  PositionCollection arbitraryPositions = 4. * Eigen::MatrixX3d::Random(arbitraryElements.size(), 3);
  auto structure = Utils::AtomCollection(arbitraryElements, arbitraryPositions);
  auto masses = Utils::Geometry::Properties::getMasses(arbitraryElements);

  // Default i.e. initial velocities drawn from 300 K Maxwell-Boltzmann distribution
  md.performMDSimulation(structure, log);
  ASSERT_THAT(md.getTemperatures()[0], Gt(50));

  // Zero velocities
  md.settings().modifyDouble(SettingsNames::generationTemperature, 0.);
  md.performMDSimulation(structure, log);
  ASSERT_THAT(md.getVelocities().size(), Eq(2));
  auto firstVelocities = md.getVelocities()[0];
  ASSERT_THAT(md.getTemperatures().size(), Eq(2));
  double firstTemperature = md.getTemperatures()[0];
  ASSERT_TRUE(firstVelocities.isZero(1e-9));
  EXPECT_THAT(firstTemperature, DoubleNear(0, 1e-9));

  // Specified velocities
  Eigen::MatrixXd randomVelocities = Eigen::MatrixXd::Random(structure.size(), 3);
  md.setInitialVelocities(randomVelocities);
  md.performMDSimulation(structure, log);
  ASSERT_THAT(md.getVelocities().size(), Eq(2));
  auto newFirstVelocities = md.getVelocities()[0];
  ASSERT_TRUE(newFirstVelocities.isApprox(randomVelocities, 1e-9));
}

TEST_F(AMolecularDynamicsTest, mdSimulationIsTerminatedEarly) {
  MockCalculator mockCalculator;
  MolecularDynamics md(mockCalculator);

  md.settings().modifyInt(SettingsNames::numberOfMDSteps, 15);
  md.settings().modifyString(SettingsNames::integrationAlgorithm, OptionNames::leapFrogOption);
  md.settings().modifyBool(SettingsNames::requireCharges, true);
  md.settings().modifyBool(SettingsNames::requireBondOrders, true);

  std::stringstream ss("3\n\n"
                       "O      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.4005871485    0.3100978455    0.0000000000\n"
                       "H     -0.4005871485    0.3100978455    0.0000000000\n");
  auto structure = Utils::XyzStreamHandler::read(ss);
  bool hasBondOrders = false;
  bool hasAtomicCharges = false;
  bool positionsCorrectSize = false;

  // Set external stop function
  auto stopFunction = [&hasBondOrders, &hasAtomicCharges, &positionsCorrectSize,
                       &structure](const PositionCollection& p, const Results& results, int steps) {
    positionsCorrectSize = p.rows() == structure.size() && p.cols() == 3;
    hasBondOrders = results.has<Property::BondOrderMatrix>();
    hasAtomicCharges = results.has<Property::AtomicCharges>();
    return steps >= 5;
  };
  md.setExternalStop(stopFunction);

  // Normal MD test
  md.performMDSimulation(structure, log);
  MolecularTrajectory mt = md.getMolecularTrajectory();

  ASSERT_THAT(mt.size(), Eq(6));
  for (const auto& energy : mt.getEnergies()) {
    ASSERT_THAT(energy, Eq(42.0));
  }
  for (const auto& pos : mt) {
    ASSERT_THAT(pos.rows(), Eq(structure.size()));
  }
  ASSERT_TRUE(hasBondOrders);
  ASSERT_TRUE(hasAtomicCharges);
  ASSERT_TRUE(positionsCorrectSize);
}

TEST_F(AMolecularDynamicsTest, mdSimulationBiasPotentialEliminatesGradient) {
  MockCalculator mockCalculator;
  MolecularDynamics md(mockCalculator);

  md.settings().modifyInt(SettingsNames::numberOfMDSteps, 10);
  md.settings().modifyString(SettingsNames::integrationAlgorithm, OptionNames::leapFrogOption);
  md.settings().modifyString(SettingsNames::thermostatAlgorithm, OptionNames::noThermostatOption);
  md.settings().modifyDouble(SettingsNames::generationTemperature, 0.0);

  std::stringstream ss("3\n\n"
                       "O      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.4005871485    0.3100978455    0.0000000000\n"
                       "H     -0.4005871485    0.3100978455    0.0000000000\n");
  auto structure = Utils::XyzStreamHandler::read(ss);
  int numSteps = 0;

  // Set bias potential
  auto bias = [&numSteps](const PositionCollection& p, const Results& results, int steps) {
    numSteps = steps;
    Utils::GradientCollection gradients(p.rows(), 3);
    gradients = -results.get<Property::Gradients>(); // the resulting gradients will be zero
    std::pair<double, GradientCollection> returnPair = std::make_pair(0.0, gradients);
    return returnPair;
  };
  md.setBiasPotential(bias);

  // Normal MD test
  md.performMDSimulation(structure, log);
  MolecularTrajectory mt = md.getMolecularTrajectory();

  ASSERT_THAT(mt.size(), Eq(11));
  for (const auto& energy : mt.getEnergies()) {
    ASSERT_THAT(energy, Eq(42.0));
  }
  for (const auto& pos : mt) {
    ASSERT_THAT(pos.rows(), Eq(structure.size()));
    ASSERT_TRUE(pos.isApprox(structure.getPositions()));
  }
  ASSERT_THAT(numSteps, Eq(10));
}

TEST_F(AMolecularDynamicsTest, mdSimulationRunsWithCalculatorWithReference) {
  auto mockCalculator = std::make_shared<MockCalculator>();
  MockCalculatorWithReference mockCalculatorWithReference;
  mockCalculatorWithReference.setReferenceCalculator(mockCalculator);
  MolecularDynamics md(mockCalculatorWithReference);

  md.settings().modifyInt(SettingsNames::numberOfMDSteps, 5);
  md.settings().modifyBool(SettingsNames::requireCharges, true);
  md.settings().modifyBool(SettingsNames::requireBondOrders, true);

  std::stringstream ss("3\n\n"
                       "O      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.4005871485    0.3100978455    0.0000000000\n"
                       "H     -0.4005871485    0.3100978455    0.0000000000\n");

  auto structure = Utils::XyzStreamHandler::read(ss);
  bool positionsCorrectSize = false;
  std::vector<double> atomicCharges;
  Utils::BondOrderCollection bondOrders;

  // Set external stop function
  // Also use this to check results contents
  auto stopFunction = [&positionsCorrectSize, &atomicCharges, &bondOrders,
                       &structure](const PositionCollection& p, const Results& results, int steps) {
    positionsCorrectSize = p.rows() == structure.size() && p.cols() == 3;
    atomicCharges = results.get<Property::AtomicCharges>();
    bondOrders = results.get<Property::BondOrderMatrix>();
    return steps >= 3;
  };
  md.setExternalStop(stopFunction);

  // Run
  md.performMDSimulation(structure, log);
  MolecularTrajectory mt = md.getMolecularTrajectory();

  ASSERT_THAT(mt.size(), Eq(4));
  for (const auto& energy : mt.getEnergies()) {
    ASSERT_THAT(energy, Eq(84.0));
  }
  for (const auto& pos : mt) {
    ASSERT_THAT(pos.rows(), Eq(structure.size()));
  }

  ASSERT_THAT(bondOrders.getSystemSize(), Eq(structure.size()));
  ASSERT_THAT(bondOrders.getOrder(0, 1), Eq(3.0));
  ASSERT_THAT(atomicCharges.size(), Eq(structure.size()));
  for (const auto& charge : atomicCharges) {
    ASSERT_THAT(charge, Eq(21.0));
  }
  ASSERT_TRUE(positionsCorrectSize);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
