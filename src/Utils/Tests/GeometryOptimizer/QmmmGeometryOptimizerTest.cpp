/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/CalculatorBasics/TestCalculator.h>
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
    Utils::UniversalSettings::IntListDescriptor qmAtoms("A list of the indices of the atoms in the QM region.");
    this->_fields.push_back("qm_atoms", std::move(qmAtoms));

    resetToDefaults();
  };
};

// Define a mock calculator
class QmmmGeoOptMockCalculator : public Core::Calculator {
 public:
  QmmmGeoOptMockCalculator() {
    this->settings_ = std::make_unique<QmmmGeoOptMockCalculatorSettings>();
    this->underlyingCalculator_ = std::make_unique<TestCalculator>();
  };
  ~QmmmGeoOptMockCalculator() override = default;
  void setStructure(const AtomCollection& structure) final {
    underlyingCalculator_->setStructure(structure);
    structure_ = structure;
  }
  void modifyPositions(PositionCollection newPositions) final {
    underlyingCalculator_->modifyPositions(newPositions);
    structure_.setPositions(newPositions);
  }
  const PositionCollection& getPositions() const final {
    return underlyingCalculator_->getPositions();
  }
  void setRequiredProperties(const PropertyList& /* requiredProperties */) final{};
  PropertyList getRequiredProperties() const final {
    return PropertyList{};
  }
  PropertyList possibleProperties() const final {
    return Utils::Property::Energy | Utils::Property::Gradients;
  }
  const Results& calculate(std::string dummy = "") final {
    results_ = underlyingCalculator_->calculate(dummy);
    return results_;
  };
  std::string name() const final {
    return std::string("QmmmGeoOptMockCalculator");
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
  std::unique_ptr<Core::Calculator> underlyingCalculator_;
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
                       "C -2.881039   -1.131901   -0.039580\n"
                       "C -3.368457   0.3026080   -0.017690\n"
                       "C -2.338492   1.4097370   0.0706880\n"
                       "C -0.855344   1.0882420   0.0874170\n"
                       "C -0.343569   -0.421243   0.0118950\n"
                       "C -1.388626   -1.412340   -0.021191\n"
                       "N 1.1791810   -0.584422   -0.064237\n"
                       "C 2.4758990   0.1334760   -0.042368\n"
                       "N 3.7426770   -0.566886   0.1428610\n"
                       "O 2.4830730   1.5113290   -0.194210\n"
                       "H 1.3528760   -1.593323   0.0005230\n"
                       "H 4.7007330   -0.227179   0.1856740\n"
                       "H 3.8773790   -1.546954   0.3937500\n"
                       "H -0.134336   1.8579500   0.0895200\n"
                       "H -1.067972   -2.422238   -0.071626\n"
                       "H -2.579605   2.4431120   0.0803140\n"
                       "H -4.393787   0.5758600   -0.006209\n"
                       "H -3.591501   -1.918491   -0.085768\n");

  // This structure is chosen such that it is not too far from the
  // minimum (w.r.t. the test calcutor).
  auto structure = Utils::XyzStreamHandler::read(ss);

  QmmmGeoOptMockCalculator calculator;
  calculator.settings().modifyIntList("qm_atoms", std::vector<int>{{8, 11, 12}});

  Core::Log log = Core::Log::silent();
  QmmmGeometryOptimizer<Bfgs> qmmmGeometryOptimizer(calculator);

  // Optimizer settings
  const int maxFullOptMicroiter = 30;
  Settings settings = qmmmGeometryOptimizer.getSettings();
  settings.modifyString(GeometryOptimizerBase::geooptCoordinateSystemKey, "cartesian");
  settings.modifyInt(QmmmGeometryOptimizer<Bfgs>::qmmmOptMaxEnvMicroiterationsKey, 100);
  settings.modifyInt(QmmmGeometryOptimizer<Bfgs>::qmmmOptMaxFullMicroiterationsKey, maxFullOptMicroiter);
  // Make sure the number of max. macroiterations is not limiting the convergence:
  settings.modifyInt(QmmmGeometryOptimizer<Bfgs>::qmmmOptMaxMacroiterationsKey, 200);
  qmmmGeometryOptimizer.setSettings(settings);

  AtomCollection s1 = structure;
  auto cyclesOne = qmmmGeometryOptimizer.optimize(s1, log);
  // Should not converge  within first round of full opt.:
  EXPECT_TRUE(cyclesOne > maxFullOptMicroiter);

  // Reset
  numberOfCalculations = 0;

  // Modified optimizer settings
  settings = qmmmGeometryOptimizer.getSettings();
  settings.modifyInt(QmmmGeometryOptimizer<Bfgs>::qmmmOptMaxEnvMicroiterationsKey, 1000);
  qmmmGeometryOptimizer.setSettings(settings);

  AtomCollection s2 = structure;
  // This optimization should now work better, because a lot more environment
  // microiterations are performed.
  auto cyclesTwo = qmmmGeometryOptimizer.optimize(s2, log);
  EXPECT_TRUE(cyclesTwo > 1); // Should converge in more than 1 cycle.
  EXPECT_TRUE(cyclesOne > cyclesTwo);

  // Reset
  numberOfCalculations = 0;

  // Modified optimizer settings
  settings = qmmmGeometryOptimizer.getSettings();
  settings.modifyInt(QmmmGeometryOptimizer<Bfgs>::qmmmOptMaxEnvMicroiterationsKey, 5);
  qmmmGeometryOptimizer.setSettings(settings);

  AtomCollection s3 = structure;
  // This optimization should do the worst of the three examples:
  // It performs only 5 environment relaxation steps in between the full opts.
  auto cyclesThree = qmmmGeometryOptimizer.optimize(s3, log);

  EXPECT_TRUE(cyclesOne < cyclesThree);
  EXPECT_TRUE(cyclesTwo < cyclesThree);

  Utils::PositionCollection optimizedPositions = s3.getPositions();
  for (int i = 0; i < optimizedPositions.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      double pos = optimizedPositions(i, j) * Constants::angstrom_per_bohr;
      EXPECT_TRUE(pos < 5.0 && pos > -5.0); // reasonable values
    }
  }

  // Final test: the optimizer should also complete the optimization when using an internal coordinate system
  settings = qmmmGeometryOptimizer.getSettings();
  settings.modifyString(GeometryOptimizerBase::geooptCoordinateSystemKey, "internal");
  qmmmGeometryOptimizer.setSettings(settings);
  qmmmGeometryOptimizer.optimize(structure, log);
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
