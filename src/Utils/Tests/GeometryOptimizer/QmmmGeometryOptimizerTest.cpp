/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/CalculatorBasics/QmmmTestCalculator.h>
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

  QmmmTestCalculator calculator;
  calculator.settings().modifyIntList("qm_atoms", std::vector<int>{{8, 11, 12}});

  Core::Log log = Core::Log::silent();
  QmmmGeometryOptimizer<Bfgs> qmmmGeometryOptimizer(calculator);

  // BFGS settings
  Settings settingsBfgs = qmmmGeometryOptimizer.fullOptimizer->getSettings();
  settingsBfgs.modifyBool(SettingsNames::Optimizations::Bfgs::useTrustRadius, false);
  qmmmGeometryOptimizer.fullOptimizer->setSettings(settingsBfgs);
  qmmmGeometryOptimizer.mmOptimizer->setSettings(settingsBfgs);
  // Optimizer settings
  const int maxFullOptMicroiter = 30;
  Settings settings = qmmmGeometryOptimizer.getSettings();
  settings.modifyString(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem, "cartesian");
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
  qmmmGeometryOptimizer.clearConstrainedAtoms();

  // Modified optimizer settings
  settings = qmmmGeometryOptimizer.getSettings();
  settings.modifyInt(QmmmGeometryOptimizer<Bfgs>::qmmmOptMaxEnvMicroiterationsKey, 1000);
  qmmmGeometryOptimizer.setSettings(settings);

  AtomCollection s2 = structure;
  // This optimization should now work better, because a lot more environment
  // microiterations are performed.
  auto cyclesTwo = qmmmGeometryOptimizer.optimize(s2, log);
  EXPECT_TRUE(cyclesTwo > 1); // Should converge in more than 1 cycle.

  // Reset
  qmmmGeometryOptimizer.clearConstrainedAtoms();

  // Modified optimizer settings
  settings = qmmmGeometryOptimizer.getSettings();
  settings.modifyInt(QmmmGeometryOptimizer<Bfgs>::qmmmOptMaxEnvMicroiterationsKey, 5);
  qmmmGeometryOptimizer.setSettings(settings);

  AtomCollection s3 = structure;
  // This optimization should do the worst of the three examples:
  // It performs only 5 environment relaxation steps in between the full opts.
  qmmmGeometryOptimizer.optimize(s3, log);

  Utils::PositionCollection optimizedPositions = s3.getPositions();
  for (int i = 0; i < optimizedPositions.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      double pos = optimizedPositions(i, j) * Constants::angstrom_per_bohr;
      EXPECT_TRUE(pos < 5.0 && pos > -5.0); // reasonable values
    }
  }

  // Final test: the optimizer should also complete the optimization when using an internal coordinate system
  settings = qmmmGeometryOptimizer.getSettings();
  settings.modifyString(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem, "internal");
  qmmmGeometryOptimizer.setSettings(settings);
  qmmmGeometryOptimizer.optimize(structure, log);
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
