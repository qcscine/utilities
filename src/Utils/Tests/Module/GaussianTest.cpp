/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/ExternalQC/ExternalProgram.h>
#include <Utils/ExternalQC/Gaussian/GaussianCalculator.h>
#include <Utils/ExternalQC/Gaussian/GaussianCalculatorSettings.h>
#include <Utils/ExternalQC/Gaussian/GaussianOutputParser.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/NativeFilenames.h>
#include <gmock/gmock.h>
#include <boost/dll/runtime_symbol_info.hpp>
#include <boost/filesystem.hpp>
#include <regex>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

class AGaussianTest : public Test {
 public:
  ExternalQC::GaussianCalculator calculator;
  boost::filesystem::path pathToResource;

 private:
  void SetUp() final {
    pathToResource = boost::dll::program_location().parent_path();
    pathToResource /= "Resources";
  }
};

TEST_F(AGaussianTest, SettingsAreSetCorrectly) {
  calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, 1);
  calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 2);
  calculator.settings().modifyString(Utils::SettingsNames::method, "PBEPBE");
  calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2SVP");
  calculator.settings().modifyString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, "test_1");
  calculator.settings().modifyString(Utils::ExternalQC::SettingsNames::gaussianFilenameBase, "test_2");
  calculator.settings().modifyInt(Utils::SettingsNames::externalProgramMemory, 4096);

  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::molecularCharge), Eq(1));
  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::spinMultiplicity), Eq(2));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::method), Eq("PBEPBE"));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::basisSet), Eq("def2SVP"));
  ASSERT_THAT(calculator.settings().getString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory), Eq("test_1"));
  ASSERT_THAT(calculator.settings().getString(Utils::ExternalQC::SettingsNames::gaussianFilenameBase), Eq("test_2"));
  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::externalProgramMemory), Eq(4096));
}

TEST_F(AGaussianTest, WorkingDirectoryCanBeSet) {
  std::string workingDir = "/home/path/to/testDir";
  ExternalQC::ExternalProgram externalProgram;
  externalProgram.setWorkingDirectory(workingDir);
  ASSERT_THAT(externalProgram.getWorkingDirectory(), NativeFilenames::addTrailingSeparator(workingDir));
}

TEST_F(AGaussianTest, OutputIsParsedCorrectly) {
  boost::filesystem::path resPath = pathToResource / "gaussian_test_calc.out";
  ExternalQC::GaussianOutputParser parser(resPath.string());

  ASSERT_THAT(parser.getEnergy(), DoubleNear(-76.27238594, 1e-4));

  Eigen::MatrixXd refGrad(3, 3);
  refGrad << -0.008031895, -0.0, -0.005389675, -0.000000000, 0.0, 0.010779349, 0.008031895, -0.0, -0.005389675;
  GradientCollection grad = parser.getGradients();

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ASSERT_THAT(grad(i, j), DoubleNear(refGrad(i, j), 1e-4));
    }
  }

  std::vector<double> refCM5Charges = {0.328269, -0.656534, 0.328269};
  std::vector<double> parsedCM5Charges = parser.getCM5Charges();

  for (int k = 0; k < refCM5Charges.size(); ++k) {
    ASSERT_THAT(parsedCM5Charges[k], DoubleNear(refCM5Charges[k], 1e-4));
  }
}

TEST_F(AGaussianTest, CloneInterfaceWorksCorrectly) {
  calculator.settings().modifyString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, pathToResource.string());

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
  ASSERT_THAT(calculator.settings().getString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory),
              Eq(newCalculator->settings().getString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory)));
}

TEST_F(AGaussianTest, InputCreationWorkCorrectly) {
  // Set up.
  calculator.settings().modifyString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, pathToResource.string());
  calculator.settings().modifyString(Utils::SettingsNames::method, "PBEPBE");
  calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2SVP FORCE_FAILURE");
  calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::AtomicCharges);

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

  std::string inputFileName =
      NativeFilenames::combinePathSegments(calculator.getCalculationDirectory(), calculator.getFileNameBase() + ".inp");

  // Check that the input file was correctly written.
  std::string line;
  std::ifstream input;
  input.open(inputFileName);
  if (input.is_open()) {
    // Get third line
    getline(input, line);
    getline(input, line);
    getline(input, line);
  }
  input.close();

  ASSERT_THAT(line, Eq("# PBEPBE/def2SVP FORCE_FAILURE Force Pop=Hirshfeld"));

  // Check whether the calculation directory can be deleted.
  bool isDir = FilesystemHelpers::isDirectory(calculator.getCalculationDirectory());
  boost::filesystem::remove_all(calculator.getCalculationDirectory());
  bool deleted = !FilesystemHelpers::isDirectory(calculator.getCalculationDirectory());
  ASSERT_THAT(isDir, Eq(true));
  ASSERT_THAT(deleted, Eq(true));
}

TEST_F(AGaussianTest, SolventInputKeywordWorkCorrectly) {
  calculator.settings().modifyString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, pathToResource.string());
  calculator.settings().modifyString(Utils::SettingsNames::solvent, "water");
  calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2SVP FORCE_FAILURE");

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

  std::string inputFileName =
      NativeFilenames::combinePathSegments(calculator.getCalculationDirectory(), calculator.getFileNameBase() + ".inp");

  std::ifstream input;
  input.open(inputFileName);
  auto content = std::string(std::istreambuf_iterator<char>{input}, {});
  input.close();

  std::string regexString = "SCRF\\=\\(CPCM\\,Solvent\\=water\\)";
  std::regex regex(regexString);
  std::smatch matches;
  bool b = std::regex_search(content, matches, regex);

  ASSERT_TRUE(b);

  // Check whether the calculation directory can be deleted.
  bool isDir = FilesystemHelpers::isDirectory(calculator.getCalculationDirectory());
  boost::filesystem::remove_all(calculator.getCalculationDirectory());
  bool deleted = !FilesystemHelpers::isDirectory(calculator.getCalculationDirectory());
  ASSERT_THAT(isDir, Eq(true));
  ASSERT_THAT(deleted, Eq(true));
}

TEST_F(AGaussianTest, GaussianCalculationIsPerformedCorrectlyViaScine) {
#ifndef _WIN32
  const char* envVariablePtr = std::getenv("GAUSSIAN_BINARY_PATH");
  if (envVariablePtr) {
    // Set-up
    calculator.settings().modifyInt(Utils::SettingsNames::externalProgramNProcs, 1);
    calculator.settings().modifyString(Utils::SettingsNames::method, "PBEPBE");
    calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2SVP");
    calculator.setRequiredProperties(Property::Energy);

    std::stringstream stream("2\n\n"
                             "H     0.00000000   0.00000000   0.00000000\n"
                             "H     0.00000000   0.00000000   0.74000000\n");
    auto structure = Utils::XyzStreamHandler::read(stream);
    calculator.setStructure(structure);

    // Calculate
    const auto& results = calculator.calculate("");
    // Assert energy
    ASSERT_THAT(results.get<Property::Energy>(), DoubleNear(-1.160019597, 1e-6));

    // Check whether the calculation directory can be deleted.
    bool isDir = FilesystemHelpers::isDirectory(calculator.getCalculationDirectory());
    boost::filesystem::remove_all(calculator.getCalculationDirectory());
    bool deleted = !FilesystemHelpers::isDirectory(calculator.getCalculationDirectory());
    ASSERT_THAT(isDir, Eq(true));
    ASSERT_THAT(deleted, Eq(true));
  }
  else {
    auto logger = Core::Log();
    logger.output << "Gaussian calculations were not tested directly as no binary path was specified." << Core::Log::endl;
  }
#endif
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
