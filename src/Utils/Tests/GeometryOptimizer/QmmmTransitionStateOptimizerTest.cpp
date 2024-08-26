/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/CalculatorBasics/QmmmEmbeddingTestCalculator.h>
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/GeometryOptimization/QmmmTransitionStateOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Settings.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class AQmmmTransitionStateOptimizerTest QmmmTransitionStateOptimizerTest.cpp
 * @brief Tests the QM/MM transition state optimizer.
 * @test
 */
class AQmmmTransitionStateOptimizerTest : public Test {
 public:
  void prepareSettingsOfUnderlyingOptimizer(std::shared_ptr<QmmmTransitionStateOptimizer<EigenvectorFollowing, Bfgs>> qmmmOptimizer) {
    // BFGS settings for MM Optimizer
    Settings settingsBfgs = qmmmOptimizer->mmOptimizer->getSettings();
    settingsBfgs.modifyBool(SettingsNames::Optimizations::Bfgs::useTrustRadius, true);
    settingsBfgs.modifyDouble(SettingsNames::Optimizations::Bfgs::trustRadius, 0.2);
    settingsBfgs.modifyDouble(SettingsNames::Optimizations::Convergence::deltaValue, 1e-2);
    settingsBfgs.modifyInt(SettingsNames::Optimizations::Convergence::requirement, 2);
    settingsBfgs.modifyDouble(SettingsNames::Optimizations::Convergence::gradMaxCoeff, 5e-3);
    settingsBfgs.modifyDouble(SettingsNames::Optimizations::Convergence::gradRMS, 1e-3);
    qmmmOptimizer->mmOptimizer->setSettings(settingsBfgs);
    // Eigenvector Following settings for QM optimizer
    Settings settingsEV = qmmmOptimizer->qmOptimizer->getSettings();
    settingsEV.modifyDouble(SettingsNames::Optimizations::EigenvectorFollowing::trustRadius, 0.05);
    settingsEV.modifyDouble(SettingsNames::Optimizations::Convergence::deltaValue, 1e-2);
    settingsEV.modifyInt(SettingsNames::Optimizations::Convergence::requirement, 1);
    settingsEV.modifyDouble(SettingsNames::Optimizations::Convergence::gradMaxCoeff, 5e-3);
    settingsEV.modifyDouble(SettingsNames::Optimizations::Convergence::gradRMS, 5e-2);
    qmmmOptimizer->qmOptimizer->setSettings(settingsEV);
  }
};

TEST_F(AQmmmTransitionStateOptimizerTest, QmmmTransitionStateOptimizerConvergenceWorksCorrectly) {
  std::stringstream ss("8\n\n"
                       "C     -0.1298977757   -0.1593511096   -0.2279748018\n"
                       "O     -1.0624371881    0.3285186148   -0.7536820217\n"
                       "H      1.1301020502    0.0307468842   -0.0756375161\n"
                       "H      0.8218997013    0.7139804431   -0.4623819386\n"
                       "H     -3.3107256895    0.3893728084   -1.4845040241\n"
                       "H     -4.0100345111    0.2407325089   -1.7560941830\n"
                       "H     -1.1660900381    2.5757574128   -2.5551010622\n"
                       "H     -1.6404768222    2.1394610516   -2.1429885433\n");

  auto structure = Utils::XyzStreamHandler::read(ss);
  std::shared_ptr<Utils::QmmmEmbeddingTestCalculator> testCalculator =
      std::make_shared<Utils::QmmmEmbeddingTestCalculator>();
  std::shared_ptr<Utils::TestCalculator> qmCalculator = std::make_shared<Utils::TestCalculator>();
  qmCalculator->setPrecision(7);
  std::shared_ptr<Utils::TestCalculator> mmCalculator = std::make_shared<Utils::TestCalculator>();
  mmCalculator->setPrecision(7);
  testCalculator->setUnderlyingCalculators(std::vector<std::shared_ptr<Core::Calculator>>{qmCalculator, mmCalculator});
  testCalculator->settings().modifyIntList("qm_atoms", std::vector<int>{{0, 1, 2, 3}});
  testCalculator->setStructure(structure);
  Core::Log log = Core::Log::silent();

  auto calculator = std::reinterpret_pointer_cast<Core::Calculator>(testCalculator);
  auto qmmmTransitionStateOptimizer =
      std::make_shared<Utils::QmmmTransitionStateOptimizer<Utils::EigenvectorFollowing, Utils::Bfgs>>(calculator);
  // QM/MM Optimizer settings
  Settings settings = qmmmTransitionStateOptimizer->getSettings();
  settings.modifyString(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem, "cartesianWithoutRotTrans");
  settings.modifyInt(qmmmTransitionStateOptimizer->qmmmOptMaxEnvMicroiterationsKey, 1000);
  settings.modifyInt(qmmmTransitionStateOptimizer->qmmmOptMaxQmMicroiterationsKey, 14);
  // Make sure the number of max. macroiterations is not limiting the convergence:
  // NOTE: current example takes 9 steps to 'convergence'
  int maxOptMacroiter = 3;
  settings.modifyInt(qmmmTransitionStateOptimizer->qmmmOptMaxMacroiterationsKey, maxOptMacroiter);
  qmmmTransitionStateOptimizer->setSettings(settings); // Overwrites the settings of the underlying optimizer
  this->prepareSettingsOfUnderlyingOptimizer(qmmmTransitionStateOptimizer);

  AtomCollection s1 = structure;
  int cyclesNotConverged = qmmmTransitionStateOptimizer->optimize(s1, log);
  // Should not converge
  EXPECT_TRUE(cyclesNotConverged == maxOptMacroiter);

  // Increase macroiterations
  maxOptMacroiter = 10;
  settings.modifyInt(qmmmTransitionStateOptimizer->qmmmOptMaxMacroiterationsKey, maxOptMacroiter);
  qmmmTransitionStateOptimizer->setSettings(settings);
  this->prepareSettingsOfUnderlyingOptimizer(qmmmTransitionStateOptimizer);
  AtomCollection s2 = structure;
  int cyclesConverged = qmmmTransitionStateOptimizer->optimize(s2, log);
  // Converges in about 90 QM-only steps
  EXPECT_TRUE(cyclesConverged < maxOptMacroiter);

  // Test cartesian coordinates
  settings.modifyString(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem, "cartesian");
  qmmmTransitionStateOptimizer->setSettings(settings);
  this->prepareSettingsOfUnderlyingOptimizer(qmmmTransitionStateOptimizer);
  AtomCollection s3 = structure;
  int cyclesCartesianConverged = qmmmTransitionStateOptimizer->optimize(s3, log);
  // Converges in about 35 QM-only steps
  EXPECT_TRUE(cyclesCartesianConverged < maxOptMacroiter);
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
