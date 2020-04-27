/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "externalQC_output_location.h"
#include <Utils/Bonds/BondDetector.h>
#include <Utils/ExternalQC/Exceptions.h>
#include <Utils/ExternalQC/ExternalProgram.h>
#include <Utils/ExternalQC/Orca/OrcaCalculator.h>
#include <Utils/ExternalQC/Orca/OrcaCalculatorSettings.h>
#include <Utils/ExternalQC/Orca/OrcaHessianOutputParser.h>
#include <Utils/ExternalQC/Orca/OrcaMainOutputParser.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/NativeFilenames.h>
#include <gmock/gmock.h>
#include <boost/filesystem.hpp>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

class AOrcaTest : public Test {
 public:
  ExternalQC::OrcaCalculator calculator;
};

TEST_F(AOrcaTest, SettingsAreSetCorrectly) {
  calculator.settings().modifyInt(Utils::ExternalQC::SettingsNames::orcaNumProcs, 2);
  calculator.settings().modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 0.0001);
  calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, 1);
  calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 2);
  calculator.settings().modifyInt(Utils::SettingsNames::maxIterations, 125);
  calculator.settings().modifyString(Utils::SettingsNames::method, "PBE D3BJ");
  calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2-SVP");
  calculator.settings().modifyString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, "test_1");
  calculator.settings().modifyString(Utils::ExternalQC::SettingsNames::orcaFilenameBase, "test_2");
  calculator.settings().modifyInt(Utils::ExternalQC::SettingsNames::externalQCMemory, 4096);

  ASSERT_THAT(calculator.settings().getInt(Utils::ExternalQC::SettingsNames::orcaNumProcs), Eq(2));
  ASSERT_THAT(calculator.settings().getDouble(Utils::SettingsNames::selfConsistanceCriterion), Eq(0.0001));
  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::molecularCharge), Eq(1));
  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::spinMultiplicity), Eq(2));
  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::maxIterations), Eq(125));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::method), Eq("PBE D3BJ"));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::basisSet), Eq("def2-SVP"));
  ASSERT_THAT(calculator.settings().getString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory), Eq("test_1"));
  ASSERT_THAT(calculator.settings().getString(Utils::ExternalQC::SettingsNames::orcaFilenameBase), Eq("test_2"));
  ASSERT_THAT(calculator.settings().getInt(Utils::ExternalQC::SettingsNames::externalQCMemory), Eq(4096));
}

TEST_F(AOrcaTest, WorkingDirectoryCanBeSet) {
  std::string workingDir = "/home/path/to/testDir";
  ExternalQC::ExternalProgram externalProgram;
  externalProgram.setWorkingDirectory(workingDir);
  ASSERT_THAT(externalProgram.getWorkingDirectory(), NativeFilenames::addTrailingSeparator(workingDir));
}

TEST_F(AOrcaTest, HessianIsParsedCorrectly) {
  ExternalQC::OrcaMainOutputParser parser(orca_test_output);
  parser.checkForErrors();
  ExternalQC::OrcaHessianOutputParser hessianParser(orca_test_hessian_output);
  HessianMatrix hessian = hessianParser.getHessian();

  ASSERT_THAT(hessian(0, 0), DoubleNear(0.36315, 1e-4));
  ASSERT_THAT(hessian(2, 3), DoubleNear(-0.28229, 1e-4));
  ASSERT_THAT(hessian(8, 8), DoubleNear(0.21312, 1e-4));
  ASSERT_THAT(hessian(5, 6), DoubleNear(0.23586, 1e-4));
  ASSERT_THAT(hessian(1, 2), DoubleNear(hessian(2, 1), 1e-6));
  ASSERT_THAT(hessian(3, 7), DoubleNear(hessian(7, 3), 1e-6));
}

TEST_F(AOrcaTest, ErrorsAreFound) {
  bool errorFound = false;
  ExternalQC::OrcaMainOutputParser parser(orca_test_error_output);
  try {
    parser.checkForErrors();
  }
  catch (ExternalQC::OutputFileParsingError) {
    errorFound = true;
  }
  ASSERT_THAT(errorFound, true);
}

TEST_F(AOrcaTest, BondOrdersAreParsedCorrectly) {
  ExternalQC::OrcaMainOutputParser parser(orca_test_bond_orders_output);
  Utils::BondOrderCollection bondOrders = parser.getBondOrders();
  double energy = parser.getEnergy();

  // Assert energy
  ASSERT_THAT(energy, DoubleNear(-154.7253016, 1e-5));
  // Some explicit bond order assertions
  ASSERT_THAT(bondOrders.getOrder(0, 1), DoubleNear(0.9582, 1e-3));
  ASSERT_THAT(bondOrders.getOrder(4, 5), DoubleNear(0.9396, 1e-3));
  ASSERT_THAT(bondOrders.getOrder(7, 8), DoubleNear(0.9915, 1e-3));

  // For this easy molecule, the bonds should be the same as evaluated by the Utils::BondDetector
  const auto& bondOrderMatrix = bondOrders.getMatrix();
  std::stringstream stream("9\n\n"
                           "H      1.8853     -0.0401      1.0854\n"
                           "C      1.2699     -0.0477      0.1772\n"
                           "H      1.5840      0.8007     -0.4449\n"
                           "H      1.5089     -0.9636     -0.3791\n"
                           "C     -0.2033      0.0282      0.5345\n"
                           "H     -0.4993     -0.8287      1.1714\n"
                           "H     -0.4235      0.9513      1.1064\n"
                           "O     -0.9394      0.0157     -0.6674\n"
                           "H     -1.8540      0.0626     -0.4252\n");

  auto structure = Utils::XyzStreamHandler::read(stream);

  // Assert size of bond order matrix
  ASSERT_THAT(bondOrders.getMatrix().cols(), Eq(structure.size()));
  ASSERT_THAT(bondOrders.getMatrix().rows(), Eq(structure.size()));

  auto otherBondOrders = Utils::BondDetector::detectBonds(structure);
  const auto& otherBondOrderMatrix = otherBondOrders.getMatrix();

  for (int i = 0; i < structure.size(); ++i) {
    for (int j = 0; j < structure.size(); ++j) {
      if (otherBondOrderMatrix.coeff(i, j) > 0.5)
        ASSERT_TRUE(bondOrderMatrix.coeff(i, j) > 0.5);
      else
        ASSERT_FALSE(bondOrderMatrix.coeff(i, j) > 0.5);
    }
  }
}

TEST_F(AOrcaTest, BondOrdersOfLargeMoleculeAreParsedCorrectly) {
  ExternalQC::OrcaMainOutputParser parser(orca_test_large_bond_orders_output);
  Utils::BondOrderCollection bondOrders = parser.getBondOrders();
  double energy = parser.getEnergy();

  // Assert energy
  ASSERT_THAT(energy, DoubleNear(-4459.355605212757, 1e-5));
  // Some explicit bond order assertions
  ASSERT_THAT(bondOrders.getOrder(0, 1), DoubleNear(0.9296, 1e-4));
  ASSERT_THAT(bondOrders.getOrder(6, 12), DoubleNear(1.2919, 1e-4));
  ASSERT_THAT(bondOrders.getOrder(161, 171), DoubleNear(0.9292, 1e-4));
  ASSERT_THAT(bondOrders.getOrder(43, 109), DoubleNear(0.8668, 1e-4));
  // Check that bondOrders of last block where parsed
  ASSERT_THAT(bondOrders.getOrder(0, 2), DoubleNear(0.6212, 1e-4));
  // Ensure that entire BondOrderMatrix was parsed
  const auto& bondOrderMatrix = bondOrders.getMatrix();
  int nAtoms = 174;
  // Assert size of bond order matrix
  ASSERT_THAT(bondOrders.getMatrix().cols(), Eq(nAtoms));
  ASSERT_THAT(bondOrders.getMatrix().rows(), Eq(nAtoms));
}

TEST_F(AOrcaTest, OutputIsParsedCorrectly) {
  ExternalQC::OrcaMainOutputParser parser(orca_test_output);

  ASSERT_THAT(parser.getEnergy(), DoubleNear(-75.818269, 1e-5));

  Eigen::MatrixXd refGrad(3, 3);
  refGrad << -0.046454155, -0.0, -0.024167757, -0.000000000, 0.0, 0.048335514, 0.046454155, -0.0, -0.024167757;
  GradientCollection grad = parser.getGradients();

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ASSERT_THAT(grad(i, j), DoubleNear(refGrad(i, j), 1e-4));
    }
  }

  // Set up thermochemistry container
  ThermochemicalContainer container;
  container.temperature = parser.getTemperature();
  container.enthalpy = parser.getEnthalpy();
  container.entropy = parser.getEntropy();
  container.zeroPointVibrationalEnergy = parser.getZeroPointVibrationalEnergy();
  container.gibbsFreeEnergy = parser.getGibbsFreeEnergy();
  container.heatCapacityP = std::numeric_limits<double>::quiet_NaN();
  container.heatCapacityV = std::numeric_limits<double>::quiet_NaN();
  ThermochemicalComponentsContainer thermochemistry;
  thermochemistry.overall = container;

  // Test thermochemistry parsing
  ASSERT_THAT(thermochemistry.overall.temperature, DoubleNear(298.15, 1e-2));
  ASSERT_THAT(thermochemistry.overall.enthalpy, DoubleNear(0.00094421 + 0.02385382, 1e-6));
  // Alternative way to calculate the enthalpy from numbers given in the output:
  ASSERT_THAT(thermochemistry.overall.enthalpy, DoubleNear(-75.79347127 - parser.getEnergy(), 1e-6));
  ASSERT_THAT(thermochemistry.overall.entropy, DoubleNear(0.02210368 / container.temperature, 1e-6));
  ASSERT_THAT(thermochemistry.overall.zeroPointVibrationalEnergy, DoubleNear(0.02101209, 1e-6));
  ASSERT_THAT(thermochemistry.overall.gibbsFreeEnergy, DoubleNear(0.00269434, 1e-6));
  ASSERT_FALSE(thermochemistry.overall.heatCapacityP == container.heatCapacityP);
  ASSERT_FALSE(thermochemistry.overall.heatCapacityV == container.heatCapacityV);
}

TEST_F(AOrcaTest, CloneInterfaceWorksCorrectly) {
  calculator.settings().modifyInt(Utils::ExternalQC::SettingsNames::orcaNumProcs, 2);
  calculator.settings().modifyString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, orca_test_dir);

  calculator.results().set<Property::Energy>(42.0);

  std::stringstream stream("5\n\n"
                           "C     0.00000000   0.00000001  -0.00000097\n"
                           "H     0.62612502   0.62612484   0.62613824\n"
                           "H    -0.62612503  -0.62612486   0.62613824\n"
                           "H    -0.62612481   0.62612463  -0.62613657\n"
                           "H     0.62612481  -0.62612464  -0.62613657\n");

  auto structure = Utils::XyzStreamHandler::read(stream);

  calculator.setStructure(structure);

  auto newCalculator = calculator.clone();

  ASSERT_THAT(calculator.getPositions()(3, 1), Eq(newCalculator->getPositions()(3, 1)));
  ASSERT_THAT(calculator.getPositions()(4, 2), Eq(newCalculator->getPositions()(4, 2)));
  ASSERT_THAT(calculator.results().get<Property::Energy>(), Eq(newCalculator->results().get<Property::Energy>()));
  ASSERT_THAT(calculator.settings().getInt(Utils::ExternalQC::SettingsNames::orcaNumProcs),
              Eq(newCalculator->settings().getInt(Utils::ExternalQC::SettingsNames::orcaNumProcs)));
  ASSERT_THAT(calculator.settings().getString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory),
              Eq(newCalculator->settings().getString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory)));
}

TEST_F(AOrcaTest, StatesHandlingAndInputCreationWorkCorrectly) {
  // Set up.
  calculator.settings().modifyString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, orca_test_dir);
  calculator.settings().modifyString(Utils::SettingsNames::method, "PBE D3BJ");
  calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2-SVP FORCE_FAILURE");

  std::stringstream stream("5\n\n"
                           "C     0.00000000   0.00000001  -0.00000097\n"
                           "H     0.62612502   0.62612484   0.62613824\n"
                           "H    -0.62612503  -0.62612486   0.62613824\n"
                           "H    -0.62612481   0.62612463  -0.62613657\n"
                           "H     0.62612481  -0.62612464  -0.62613657\n");

  auto structure = Utils::XyzStreamHandler::read(stream);

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

  auto state = calculator.getState();

  gbw.open(gbwFileName);
  if (gbw.is_open()) {
    gbw << "Changed Content" << std::endl;
  }
  gbw.close();

  calculator.loadState(state);

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

  ASSERT_THAT(line, Eq("! PBE D3BJ def2-SVP FORCE_FAILURE"));

  // Check whether the calculation directory can be deleted.
  bool isDir = FilesystemHelpers::isDirectory(calculator.getCalculationDirectory());
  boost::filesystem::remove_all(calculator.getCalculationDirectory());
  bool deleted = !FilesystemHelpers::isDirectory(calculator.getCalculationDirectory());
  ASSERT_THAT(isDir, Eq(true));
  ASSERT_THAT(deleted, Eq(true));
}

TEST_F(AOrcaTest, TemporaryFilesAreDeletedCorrectly) {
  // Set up.
  calculator.settings().modifyString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, orca_test_dir);
  calculator.settings().modifyString(Utils::SettingsNames::method, "FORCE_FAILURE");

  std::stringstream stream("5\n\n"
                           "C     0.00000000   0.00000001  -0.00000097\n"
                           "H     0.62612502   0.62612484   0.62613824\n"
                           "H    -0.62612503  -0.62612486   0.62613824\n"
                           "H    -0.62612481   0.62612463  -0.62613657\n"
                           "H     0.62612481  -0.62612464  -0.62613657\n");

  auto structure = Utils::XyzStreamHandler::read(stream);

  calculator.setStructure(structure);

  // Define calculation directory
  boost::filesystem::path calcDirPath(calculator.getCalculationDirectory());
  const auto& calcDir = calcDirPath.string();

  // Create calculation directory
  if (!calculator.getCalculationDirectory().empty())
    FilesystemHelpers::createDirectories(calcDir);

  // Create some temporary files
  std::string filename1 = NativeFilenames::combinePathSegments(calcDir, "1.tmp");
  std::string filename2 = NativeFilenames::combinePathSegments(calcDir, "2.tmp");
  std::string filename3 = NativeFilenames::combinePathSegments(calcDir, "abc.tmp");
  std::string filenameToKeep = NativeFilenames::combinePathSegments(calcDir, "keep_this_file.txt");
  std::ofstream file1{filename1};
  file1.close();
  std::ofstream file2{filename2};
  file2.close();
  std::ofstream file3{filename3};
  file3.close();
  std::ofstream fileToKeep{filenameToKeep};
  fileToKeep.close();

  try {
    calculator.calculate("");
  }
  catch (Core::UnsuccessfulCalculationException& e) {
  }

  // Check that all of the .tmp files were correctly deleted
  ASSERT_FALSE(boost::filesystem::exists(filename1));
  ASSERT_FALSE(boost::filesystem::exists(filename2));
  ASSERT_FALSE(boost::filesystem::exists(filename3));

  // Check that the .txt file still exists
  ASSERT_TRUE(boost::filesystem::exists(filenameToKeep));

  // Repeat the same procedure, but with the tmp file removal setting disabled
  calculator.settings().modifyBool(Utils::ExternalQC::SettingsNames::deleteTemporaryFiles, false);

  // Create the temporary files again
  std::string filename4 = NativeFilenames::combinePathSegments(calcDir, "3.tmp");
  std::string filename5 = NativeFilenames::combinePathSegments(calcDir, "4.tmp");
  std::string filename6 = NativeFilenames::combinePathSegments(calcDir, "xyz.tmp");
  std::ofstream file4{filename4};
  file4.close();
  std::ofstream file5{filename5};
  file5.close();
  std::ofstream file6{filename6};
  file6.close();

  // Check that all of the .tmp files still exist
  ASSERT_TRUE(boost::filesystem::exists(filename4));
  ASSERT_TRUE(boost::filesystem::exists(filename5));
  ASSERT_TRUE(boost::filesystem::exists(filename6));

  // Check that the .txt file still exists
  ASSERT_TRUE(boost::filesystem::exists(filenameToKeep));

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
