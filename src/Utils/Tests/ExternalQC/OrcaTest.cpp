/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "orca_output_location.h"
#include <Utils/ExternalQC/ExternalProgram.h>
#include <Utils/ExternalQC/Orca/OrcaCalculator.h>
#include <Utils/ExternalQC/Orca/OrcaCalculatorSettings.h>
#include <Utils/ExternalQC/Orca/OrcaHessianOutputParser.h>
#include <Utils/ExternalQC/Orca/OrcaMainOutputParser.h>
#include <Utils/IO/ChemicalFileFormats/XYZStreamHandler.h>
#include <Utils/IO/NativeFilenames.h>
#include <gmock/gmock.h>
#include <boost/filesystem.hpp>
#include <fstream>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

class AOrcaTest : public Test {
 public:
  ExternalQC::OrcaCalculator calculator;
};

TEST_F(AOrcaTest, SettingsAreSetCorrectly) {
  calculator.settings().modifyInt(Utils::SettingsNames::orcaNumProcs, 2);
  calculator.settings().modifyString(Utils::SettingsNames::orcaBinaryPath, "ORCAPATH");
  calculator.settings().modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 0.0001);
  calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, 1);
  calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 2);
  calculator.settings().modifyInt(Utils::SettingsNames::maxIterations, 125);
  calculator.settings().modifyString(Utils::SettingsNames::orcaMethod, "PBE D3BJ def2-SVP");
  calculator.settings().modifyString(Utils::SettingsNames::baseWorkingDirectory, "test_1");
  calculator.settings().modifyString(Utils::SettingsNames::orcaFilenameBase, "test_2");

  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::orcaNumProcs), Eq(2));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::orcaBinaryPath), Eq("ORCAPATH"));
  ASSERT_THAT(calculator.settings().getDouble(Utils::SettingsNames::selfConsistanceCriterion), Eq(0.0001));
  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::molecularCharge), Eq(1));
  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::spinMultiplicity), Eq(2));
  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::maxIterations), Eq(125));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::orcaMethod), Eq("PBE D3BJ def2-SVP"));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::baseWorkingDirectory), Eq("test_1"));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::orcaFilenameBase), Eq("test_2"));
}

TEST_F(AOrcaTest, WorkingDirectoryCanBeSet) {
  std::string workingDir = "/home/path/to/testDir";
  ExternalQC::ExternalProgram externalProgram;
  externalProgram.setWorkingDirectory(workingDir);
  ASSERT_THAT(externalProgram.getWorkingDirectory(), NativeFilenames::addTrailingSeparator(workingDir));
}

TEST_F(AOrcaTest, HessianIsParsedCorrectly) {
  ExternalQC::OrcaHessianOutputParser hessianParser(orca_test_hessian_output);
  HessianMatrix hessian = hessianParser.getHessian();

  ASSERT_THAT(hessian(0, 0), DoubleNear(0.36315, 1e-4));
  ASSERT_THAT(hessian(2, 3), DoubleNear(-0.28229, 1e-4));
  ASSERT_THAT(hessian(8, 8), DoubleNear(0.21312, 1e-4));
  ASSERT_THAT(hessian(5, 6), DoubleNear(0.23586, 1e-4));
  ASSERT_THAT(hessian(1, 2), DoubleNear(hessian(2, 1), 1e-6));
  ASSERT_THAT(hessian(3, 7), DoubleNear(hessian(7, 3), 1e-6));
}

TEST_F(AOrcaTest, OutputIsParsedCorrectly) {
  ExternalQC::OrcaMainOutputParser parser(orca_test_output);

  ASSERT_THAT(parser.getEnergy(), DoubleNear(-75.818269, 1e-4));

  Eigen::MatrixXd refGrad(3, 3);
  refGrad << -0.046454155, -0.0, -0.024167757, -0.000000000, 0.0, 0.048335514, 0.046454155, -0.0, -0.024167757;
  GradientCollection grad = parser.getGradients();

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ASSERT_THAT(grad(i, j), DoubleNear(refGrad(i, j), 1e-4));
    }
  }
}

TEST_F(AOrcaTest, CloneInterfaceWorksCorrectly) {
  calculator.settings().modifyInt(Utils::SettingsNames::orcaNumProcs, 2);
  calculator.settings().modifyString(Utils::SettingsNames::baseWorkingDirectory, orca_test_dir);

  calculator.results().setEnergy(42.0);

  std::stringstream stream("5\n\n"
                           "C     0.00000000   0.00000001  -0.00000097\n"
                           "H     0.62612502   0.62612484   0.62613824\n"
                           "H    -0.62612503  -0.62612486   0.62613824\n"
                           "H    -0.62612481   0.62612463  -0.62613657\n"
                           "H     0.62612481  -0.62612464  -0.62613657\n");

  auto structure = Utils::XYZStreamHandler::read(stream);

  calculator.setStructure(structure);

  auto newCalculator = calculator.clone();

  ASSERT_THAT(calculator.getPositions()(3, 1), Eq(newCalculator->getPositions()(3, 1)));
  ASSERT_THAT(calculator.getPositions()(4, 2), Eq(newCalculator->getPositions()(4, 2)));
  ASSERT_THAT(calculator.results().getEnergy(), Eq(newCalculator->results().getEnergy()));
  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::orcaNumProcs),
              Eq(newCalculator->settings().getInt(Utils::SettingsNames::orcaNumProcs)));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::baseWorkingDirectory),
              Eq(newCalculator->settings().getString(Utils::SettingsNames::baseWorkingDirectory)));
}

TEST_F(AOrcaTest, StatesHandlingAndInputCreationWorkCorrectly) {
  // Set up.
  calculator.settings().modifyString(Utils::SettingsNames::baseWorkingDirectory, orca_test_dir);
  calculator.settings().modifyString(Utils::SettingsNames::orcaBinaryPath, "PLACEHOLDER");
  calculator.settings().modifyString(Utils::SettingsNames::orcaMethod, "PBE D3BJ def2-SVP");

  std::stringstream stream("5\n\n"
                           "C     0.00000000   0.00000001  -0.00000097\n"
                           "H     0.62612502   0.62612484   0.62613824\n"
                           "H    -0.62612503  -0.62612486   0.62613824\n"
                           "H    -0.62612481   0.62612463  -0.62613657\n"
                           "H     0.62612481  -0.62612464  -0.62613657\n");

  auto structure = Utils::XYZStreamHandler::read(stream);

  calculator.setStructure(structure);

  try {
    calculator.calculate("");
  }
  catch (Core::UnsuccessfulCalculationException& e) {
  }

  // Check that the states handler works.
  std::string testString = "Test";
  std::string gbwFileName =
      NativeFilenames::combinePathSegments(calculator.getCalculationDirectory(), calculator.getFileNameBase() + ".gbw");
  std::string inputFileName =
      NativeFilenames::combinePathSegments(calculator.getCalculationDirectory(), calculator.getFileNameBase() + ".inp");
  std::ofstream gbw;
  gbw.open(gbwFileName);
  if (gbw.is_open()) {
    gbw << testString << std::endl;
  }
  gbw.close();

  calculator.statesHandler().store(Utils::StateSize::minimal);

  gbw.open(gbwFileName);
  if (gbw.is_open()) {
    gbw << "Changed Content" << std::endl;
  }
  gbw.close();

  calculator.statesHandler().load(calculator.statesHandler().popNewestState());

  std::string gbwContent;
  std::ifstream check;
  check.open(gbwFileName);
  if (check.is_open()) {
    check >> gbwContent;
  }
  check.close();

  ASSERT_THAT(gbwContent, Eq(testString));

  // Check that the input file was correctly written.
  std::string line;
  std::ifstream input;
  input.open(inputFileName);
  if (input.is_open()) {
    getline(input, line);
  }
  input.close();

  ASSERT_THAT(line, Eq("! PBE D3BJ def2-SVP"));

  // Check whether the calculation directory can be deleted.
  bool isDir = FilesystemHelpers::isDirectory(calculator.getCalculationDirectory());
  boost::filesystem::remove_all(calculator.getCalculationDirectory());
  bool deleted = !FilesystemHelpers::isDirectory(calculator.getCalculationDirectory());
  ASSERT_THAT(isDir, Eq(true));
  ASSERT_THAT(deleted, Eq(true));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
