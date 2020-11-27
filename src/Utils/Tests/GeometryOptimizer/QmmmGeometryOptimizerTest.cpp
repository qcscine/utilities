/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/GeometryOptimization/QmmmGeometryOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Optimizer/GradientBased/GradientBasedCheck.h>
#include <Utils/Settings.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

namespace {
// Counter for the number of calculations to adapt energy and gradients in mock calculator.
int numberOfCalculations = 0;
} // namespace

// Define mock settings
class QmmmGeoOptMockCalculatorSettings : public Scine::Utils::Settings {
 public:
  QmmmGeoOptMockCalculatorSettings() : Settings("QmmmGeoOptMockCalculatorSettings") {
    Utils::UniversalSettings::BoolDescriptor ignoreQm(
        "Whether to ignore all contributions from the QM calculation, and therefore, not performing it.");
    ignoreQm.setDefaultValue(false);
    this->_fields.push_back("ignore_qm", std::move(ignoreQm));
    Utils::UniversalSettings::StringDescriptor qmAtoms(
        "A string containing the indices of the atoms in the QM region delimited by spaces.");
    qmAtoms.setDefaultValue("");
    this->_fields.push_back("qm_atoms", std::move(qmAtoms));

    resetToDefaults();
  };
};

// Define a mock calculator
class QmmmGeoOptMockCalculator : public Core::Calculator {
 public:
  QmmmGeoOptMockCalculator() {
    this->settings_ = std::make_unique<QmmmGeoOptMockCalculatorSettings>();
  };
  ~QmmmGeoOptMockCalculator() = default;
  void setStructure(const AtomCollection& structure) final {
    structure_ = structure;
  }
  void modifyPositions(PositionCollection newPositions) final {
    structure_.setPositions(newPositions);
  }
  const PositionCollection& getPositions() const final {
    return structure_.getPositions();
  }
  void setRequiredProperties(const PropertyList& requiredProperties) final{};
  PropertyList getRequiredProperties() const final {
    return PropertyList{};
  }
  PropertyList possibleProperties() const final {
    return Utils::Property::Energy | Utils::Property::Gradients;
  }
  const Results& calculate(std::string dummy = "") final {
    numberOfCalculations++;

    results_.set<Property::Description>(dummy);

    // Harmonic potential gets less steep over time
    double forceConstant = 1e-2 / std::pow(numberOfCalculations, 2);
    if (numberOfCalculations > 100)
      forceConstant = 1e-6;

    auto energy = 0.0;
    Utils::PositionCollection positions = structure_.getPositions();
    for (int i = 0; i < positions.rows(); ++i) {
      for (int j = 0; j < 3; ++j) {
        double coord = positions(i, j);
        // Harmonic potential:
        energy += 0.5 * forceConstant * (coord - 5) * (coord - 5);
      }
    }

    results_.set<Property::Energy>(energy);
    GradientCollection grad(structure_.size(), 3);

    for (int i = 0; i < grad.rows(); ++i) {
      for (int j = 0; j < 3; ++j) {
        double coord = positions(i, j);
        grad(i, j) = -forceConstant * (coord - 5);
      }
    }

    results_.set<Property::Gradients>(grad);
    results_.set<Property::SuccessfulCalculation>(true);
    return results_;
  };
  std::string name() const final {
    return std::string("QmmmGeoOptMockCalculator");
  };
  std::shared_ptr<Core::State> getState() const final {
    return nullptr;
  }
  void loadState(std::shared_ptr<Core::State> state) final {
  }
  const Settings& settings() const final {
    return *settings_;
  }
  Settings& settings() final {
    return *settings_;
  }
  Utils::Results& results() final {
    return results_;
  }
  const Utils::Results& results() const final {
    return results_;
  }
  std::unique_ptr<Utils::AtomCollection> getStructure() const final {
    return std::make_unique<Utils::AtomCollection>(structure_);
  }
  bool supportsMethodFamily(const std::string& methodFamily) const final {
    return methodFamily == "QMMM";
  }

 private:
  AtomCollection structure_;
  Results results_;
  std::unique_ptr<Settings> settings_;
  Core::Calculator* cloneImpl() const final {
    return nullptr;
  }
};

/**
 * @class AQmmmGeometryOptimizerTest QmmmGeometryOptimizerTest.cpp
 * @brief Tests the QM/MM geometry optimizer.
 * @test
 */
class AQmmmGeometryOptimizerTest : public Test {};

TEST_F(AQmmmGeometryOptimizerTest, QmmmGeometryOptimizerConvergenceWorksCorrectly) {
  std::stringstream ss("18\n\n"
                       "C      -2.585060   -1.031121   -0.039836\n"
                       "C      -3.019151    0.286906    0.001659\n"
                       "C      -2.077282    1.304388    0.048828\n"
                       "C      -0.716691    1.033287    0.054999\n"
                       "C      -0.285145   -0.293334    0.012905\n"
                       "C      -1.232035   -1.319904   -0.034006\n"
                       "N       1.063351   -0.670447    0.030124\n"
                       "C       2.181849    0.137092   -0.019814\n"
                       "N       3.362904   -0.587503   -0.014599\n"
                       "O       2.161817    1.345506   -0.105228\n"
                       "H       1.221439   -1.662067   -0.028326\n"
                       "H       4.173941   -0.013589    0.137055\n"
                       "H       3.393655   -1.474935    0.457432\n"
                       "H       0.009797    1.828633    0.085145\n"
                       "H      -0.897370   -2.351338   -0.066181\n"
                       "H      -2.401084    2.337170    0.080642\n"
                       "H      -4.076370    0.515908   -0.003521\n"
                       "H      -3.301242   -1.842064   -0.077515\n");

  auto structure = Utils::XyzStreamHandler::read(ss);

  QmmmGeoOptMockCalculator calculator;
  calculator.settings().modifyString("qm_atoms", "8 11 12");

  Core::Log log = Core::Log::silent();
  QmmmGeometryOptimizer<Bfgs> qmmmGeometryOptimizer(calculator);

  // Optimizer settings
  const int maxFullOptMicroiter = 10;
  auto settings = qmmmGeometryOptimizer.getSettings();
  settings.modifyBool(GeometryOptimizerBase::geooptTransformCoordinatesKey, false);
  settings.modifyInt(QmmmGeometryOptimizer<Bfgs>::qmmmOptMaxEnvMicroiterationsKey, 20);
  settings.modifyInt(QmmmGeometryOptimizer<Bfgs>::qmmmOptMaxFullMicroiterationsKey, maxFullOptMicroiter);
  settings.modifyDouble(GradientBasedCheck::gconvDeltaValueKey, 1e-9); // tighter convergence criteria
  // Make sure the number of max. macroiterations is not limiting the convergence:
  settings.modifyInt(QmmmGeometryOptimizer<Bfgs>::qmmmOptMaxMacroiterationsKey, 200);
  qmmmGeometryOptimizer.setSettings(settings);

  AtomCollection s1 = structure;
  auto cyclesOne = qmmmGeometryOptimizer.optimize(s1, log);
  // Should not converge  within first round of full opt.:
  ASSERT_TRUE(cyclesOne > maxFullOptMicroiter);

  // Reset
  numberOfCalculations = 0;

  // Modified optimizer settings
  settings = qmmmGeometryOptimizer.getSettings();
  settings.modifyInt(QmmmGeometryOptimizer<Bfgs>::qmmmOptMaxEnvMicroiterationsKey, 500);
  qmmmGeometryOptimizer.setSettings(settings);

  AtomCollection s2 = structure;
  // This optimization should now work better, because a lot more environment
  // microiterations are performed.
  auto cyclesTwo = qmmmGeometryOptimizer.optimize(s2, log);
  // Should converge within first round of full opt.:
  ASSERT_TRUE(cyclesTwo < maxFullOptMicroiter);
  ASSERT_TRUE(cyclesTwo > 1); // Should converge in more than 1 cycle.
  ASSERT_TRUE(cyclesOne > cyclesTwo);

  // Reset
  numberOfCalculations = 0;

  // Modified optimizer settings
  settings = qmmmGeometryOptimizer.getSettings();
  settings.modifyInt(QmmmGeometryOptimizer<Bfgs>::qmmmOptMaxEnvMicroiterationsKey, 3);
  qmmmGeometryOptimizer.setSettings(settings);

  AtomCollection s3 = structure;
  // This optimization should do the worst of the three examples:
  // It performs only 3 environment relaxation steps in between the full opts.
  auto cyclesThree = qmmmGeometryOptimizer.optimize(s3, log);
  ASSERT_TRUE(cyclesOne < cyclesThree);
  ASSERT_TRUE(cyclesTwo < cyclesThree);

  Utils::PositionCollection optimizedPositions = s3.getPositions();
  for (int i = 0; i < optimizedPositions.rows(); ++i) {
    for (int j = 0; j < 3; ++j)
      ASSERT_THAT(optimizedPositions(i, j), DoubleNear(5.0, 1e-3));
  }
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
