/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Bonds/BondDetector.h>
#include <Utils/ExternalQC/Exceptions.h>
#include <Utils/ExternalQC/Turbomole/TurbomoleCalculator.h>
#include <Utils/ExternalQC/Turbomole/TurbomoleCalculatorSettings.h>
#include <Utils/ExternalQC/Turbomole/TurbomoleFiles.h>
#include <Utils/ExternalQC/Turbomole/TurbomoleHelper.h>
#include <Utils/ExternalQC/Turbomole/TurbomoleMainOutputParser.h>
#include <Utils/Geometry/ElementTypes.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/NativeFilenames.h>
#include <gmock/gmock.h>
#include <math.h>
#include <boost/dll/runtime_symbol_info.hpp>
#include <boost/filesystem.hpp>
#include <regex>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

class ATurbomoleTest : public Test {
 public:
  ExternalQC::TurbomoleCalculator calculator;
  ExternalQC::TurbomoleFiles files;
  std::string ressourcesDir;

 private:
  void SetUp() final {
    boost::filesystem::path pathToResource = boost::dll::program_location().parent_path();
    pathToResource /= "Resources";
    ressourcesDir = NativeFilenames::combinePathSegments(pathToResource.string(), "Turbomole");
    setCorrectTurbomoleFileNames(files, ressourcesDir);
  }
};

TEST_F(ATurbomoleTest, SettingsAreSetCorrectly) {
  calculator.settings().modifyInt(Utils::SettingsNames::externalProgramNProcs, 2);
  calculator.settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 0.0001);
  calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, 1);
  calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 2);
  calculator.settings().modifyInt(Utils::SettingsNames::maxScfIterations, 125);
  calculator.settings().modifyString(Utils::SettingsNames::method, "PBE-D3BJ");
  calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2-SVP");
  calculator.settings().modifyString(Utils::SettingsNames::spinMode, "unrestricted");
  calculator.settings().modifyString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, "test_1");
  calculator.settings().modifyDouble(Utils::SettingsNames::temperature, 300.3);
  calculator.settings().modifyString(Utils::SettingsNames::solvation, "cosmo");
  calculator.settings().modifyString(Utils::SettingsNames::solvent, "water");

  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::externalProgramNProcs), Eq(2));
  ASSERT_THAT(calculator.settings().getDouble(Utils::SettingsNames::selfConsistenceCriterion), Eq(0.0001));
  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::molecularCharge), Eq(1));
  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::spinMultiplicity), Eq(2));
  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::maxScfIterations), Eq(125));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::method), Eq("PBE-D3BJ"));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::basisSet), Eq("def2-SVP"));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::spinMode), Eq("unrestricted"));
  ASSERT_THAT(calculator.settings().getString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory), Eq("test_1"));
  ASSERT_THAT(calculator.settings().getDouble(Utils::SettingsNames::temperature), Eq(300.3));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::solvation), Eq("cosmo"));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::solvent), Eq("water"));
}

TEST_F(ATurbomoleTest, HessianOutputIsParsedCorrectly) {
  ExternalQC::TurbomoleMainOutputParser parser(files);

  HessianMatrix hessian = parser.getHessian();

  ASSERT_THAT(hessian(0, 0), DoubleNear(0.62348, 1e-4));
  ASSERT_THAT(hessian(2, 3), DoubleNear(0.0, 1e-4));
  ASSERT_THAT(hessian(7, 7), DoubleNear(0.24548, 1e-4));
  ASSERT_THAT(hessian(4, 6), DoubleNear(-0.03032, 1e-4));
  ASSERT_THAT(hessian(6, 4), DoubleNear(hessian(4, 6), 1e-6));
  ASSERT_THAT(hessian(3, 2), DoubleNear(hessian(2, 3), 1e-6));
}

TEST_F(ATurbomoleTest, NumberOfAtomsIsParsedCorrectly) {
  ExternalQC::TurbomoleMainOutputParser parser(files);
  int nAtoms = parser.getNumberAtoms();
  ASSERT_EQ(nAtoms, 3);
}

TEST_F(ATurbomoleTest, EnergyAndChargesAreParsedCorrectly) {
  ExternalQC::TurbomoleMainOutputParser parser(files);

  double energy = parser.getEnergy();
  ASSERT_THAT(energy, DoubleNear(-76.2655656, 1e-7));

  auto charges = parser.getLoewdinCharges();

  std::vector<double> correctCharges = {-0.70130, 0.35089, 0.35041};
  ASSERT_THAT(correctCharges.size(), Eq(charges.size()));
  for (unsigned i = 0; i < charges.size(); ++i) {
    ASSERT_THAT(charges[i], DoubleNear(correctCharges[i], 1e-5));
  }
}

TEST_F(ATurbomoleTest, BondOdersAreParsedCorrectly) {
  ExternalQC::TurbomoleMainOutputParser parser(files);
  Utils::BondOrderCollection bondOrders = parser.getBondOrders();

  // Some explicit bond order assertions
  ASSERT_THAT(bondOrders.getOrder(0, 1), DoubleNear(0.86939, 1e-5));
  ASSERT_THAT(bondOrders.getOrder(0, 2), DoubleNear(0.86948, 1e-5));

  // For this easy molecule, the bonds should be the same as evaluated by the Utils::BondDetector
  const auto& bondOrderMatrix = bondOrders.getMatrix();
  std::stringstream stream("3\n\n"
                           "O     0.0400775   -0.3952191    0.0000000\n"
                           "H    -0.7828149    0.1203164    0.0000000\n"
                           "H     0.7427374    0.2749027    0.0000000\n");

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

TEST_F(ATurbomoleTest, CheckResultsClearing1) {
  calculator.results().set<Property::Energy>(42.0);
  std::stringstream stream("5\n\n"
                           "C     0.00000000   0.00000001  -0.00000097\n"
                           "H     0.62612502   0.62612484   0.62613824\n"
                           "H    -0.62612503  -0.62612486   0.62613824\n"
                           "H    -0.62612481   0.62612463  -0.62613657\n"
                           "H     0.62612481  -0.62612464  -0.62613657\n");

  auto structure = Utils::XyzStreamHandler::read(stream);
  calculator.setStructure(structure);
  ASSERT_FALSE(calculator.results().has<Property::Energy>());
}

TEST_F(ATurbomoleTest, CheckResultsClearing2) {
  std::stringstream stream("5\n\n"
                           "C     0.00000000   0.00000001  -0.00000097\n"
                           "H     0.62612502   0.62612484   0.62613824\n"
                           "H    -0.62612503  -0.62612486   0.62613824\n"
                           "H    -0.62612481   0.62612463  -0.62613657\n"
                           "H     0.62612481  -0.62612464  -0.62613657\n");

  auto structure = Utils::XyzStreamHandler::read(stream);
  calculator.setStructure(structure);
  calculator.results().set<Property::Energy>(42.0);
  calculator.modifyPositions(structure.getPositions());
  ASSERT_FALSE(calculator.results().has<Property::Energy>());
}

TEST_F(ATurbomoleTest, CloneInterfaceWorksCorrectly) {
  calculator.settings().modifyInt(Utils::SettingsNames::externalProgramNProcs, 2);
  calculator.settings().modifyString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, ressourcesDir);

  std::stringstream stream("5\n\n"
                           "C     0.00000000   0.00000001  -0.00000097\n"
                           "H     0.62612502   0.62612484   0.62613824\n"
                           "H    -0.62612503  -0.62612486   0.62613824\n"
                           "H    -0.62612481   0.62612463  -0.62613657\n"
                           "H     0.62612481  -0.62612464  -0.62613657\n");

  auto structure = Utils::XyzStreamHandler::read(stream);

  calculator.setStructure(structure);

  calculator.results().set<Property::Energy>(42.0);

  auto newCalculator = calculator.clone();

  ASSERT_THAT(calculator.getPositions()(3, 1), Eq(newCalculator->getPositions()(3, 1)));
  ASSERT_THAT(calculator.getPositions()(4, 2), Eq(newCalculator->getPositions()(4, 2)));
  ASSERT_THAT(calculator.results().get<Property::Energy>(), Eq(newCalculator->results().get<Property::Energy>()));
  ASSERT_THAT(calculator.settings().getInt(Utils::SettingsNames::externalProgramNProcs),
              Eq(newCalculator->settings().getInt(Utils::SettingsNames::externalProgramNProcs)));
  ASSERT_THAT(calculator.settings().getString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory),
              Eq(newCalculator->settings().getString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory)));
}

TEST_F(ATurbomoleTest, TurbomoleCalculationIsPerformedCorrectlyViaScine) {
#ifndef _WIN32
  const char* envVariablePtr = std::getenv("TURBODIR");
  if (envVariablePtr) {
    calculator.settings().modifyInt(Utils::SettingsNames::externalProgramNProcs, 1);
    calculator.settings().modifyString(Utils::SettingsNames::method, "pbe d3bj");
    calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
    calculator.setRequiredProperties(Property::Energy);

    std::stringstream stream("5\n\n"
                             "C     0.00000000   0.00000001  -0.00000097\n"
                             "H     0.62612502   0.62612484   0.62613824\n"
                             "H    -0.62612503  -0.62612486   0.62613824\n"
                             "H    -0.62612481   0.62612463  -0.62613657\n"
                             "H     0.62612481  -0.62612464  -0.62613657\n");
    auto structure = Utils::XyzStreamHandler::read(stream);
    calculator.setStructure(structure);
    calculator.setRequiredProperties(Property::Energy);

    // dispersion correction is given incorrectly
    calculator.settings().modifyString(Utils::SettingsNames::method, "pbe D3D3");
    ASSERT_THROW(calculator.calculate(""), std::logic_error);
    calculator.settings().modifyString(Utils::SettingsNames::method, "pbe-d3bj");

    // charge and multiplicity pair does not match
    calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 2);
    calculator.settings().modifyString(Utils::SettingsNames::spinMode, "unrestricted");
    calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, 0);
    ASSERT_THROW(calculator.calculate(""), std::logic_error);
    calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 1);

    // spin mode is given incorrectly
    calculator.settings().modifyString(Utils::SettingsNames::spinMode, "restricted_open_shell");
    ASSERT_THROW(calculator.calculate(""), std::logic_error);
    calculator.settings().modifyString(Utils::SettingsNames::spinMode, "unrestricted");

    // solvent is not implemented
    calculator.settings().modifyString(Utils::SettingsNames::solvation, "cosmo");
    calculator.settings().modifyString(Utils::SettingsNames::solvent, "methane");
    ASSERT_THROW(calculator.calculate(""), std::runtime_error);
    calculator.settings().modifyString(Utils::SettingsNames::solvation, "");
    calculator.settings().modifyString(Utils::SettingsNames::solvent, "");

    // unconverged SCF is detected
    calculator.settings().modifyInt(Utils::SettingsNames::maxScfIterations, 2);
    ASSERT_THROW(calculator.calculate(""), std::runtime_error);
    calculator.settings().modifyInt(Utils::SettingsNames::maxScfIterations, 100);

    // Calculate
    const auto& results = calculator.calculate("");

    // Assert energy
    ASSERT_THAT(results.get<Property::Energy>(), DoubleNear(-40.41523066560, 1e-8));

    // Check whether the calculation directory can be deleted.
    bool isDir = FilesystemHelpers::isDirectory(calculator.getCalculationDirectory());
    boost::filesystem::remove_all(calculator.getCalculationDirectory());
    bool deleted = !FilesystemHelpers::isDirectory(calculator.getCalculationDirectory());
    ASSERT_THAT(isDir, Eq(true));
    ASSERT_THAT(deleted, Eq(true));
  }
  else {
    auto logger = Core::Log();
    logger.output << "Turbomole calculations were not tested directly as no binary path was specified." << Core::Log::endl;
  }
#endif
}

TEST_F(ATurbomoleTest, GradientsAreParsedCorrectly) {
  ExternalQC::TurbomoleMainOutputParser parser(files);
  GradientCollection grad = parser.getGradients();

  Eigen::MatrixXd refGrad(3, 3);

  refGrad << -0.00143449, 0.0140125, 1.42904e-06, 0.00660382, -0.00640583, -4.81588e-06, -0.00522462, -0.00767626, 3.387e-06;

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ASSERT_THAT(grad(i, j), DoubleNear(refGrad(i, j), 1e-6));
    }
  }

  GradientCollection pcGrad = parser.getPointChargesGradients();

  Eigen::MatrixXd refPcGrad(3, 3);
  refPcGrad << -0.000466463, 1.30447e-05, 0.000307525, -0.000362372, -5.79554e-05, 0.000517804, 0.000865896,
      -4.4095e-05, -0.000890785;

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ASSERT_THAT(pcGrad(i, j), DoubleNear(refPcGrad(i, j), 1e-8));
    }
  }
}

TEST_F(ATurbomoleTest, InputFileIsWrittenCorrectlyAndStatesHandlingWorks) {
#ifndef _WIN32
  const char* envVariablePtr = std::getenv("TURBODIR");
  if (envVariablePtr) {
    calculator.settings().modifyInt(Utils::SettingsNames::maxScfIterations, 1);
    calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, 1);
    calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 2);
    calculator.settings().modifyString(Utils::SettingsNames::solvation, "cosmo");
    calculator.settings().modifyString(Utils::SettingsNames::solvent, "water");

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

    auto calcDir = calculator.getCalculationDirectory();
    ExternalQC::TurbomoleFiles outputFiles;
    setCorrectTurbomoleFileNames(outputFiles, calcDir);

    ASSERT_TRUE(boost::filesystem::exists(outputFiles.controlFile));
    ASSERT_TRUE(boost::filesystem::exists(outputFiles.defineInputFile));
    ASSERT_TRUE(boost::filesystem::exists(outputFiles.coordFile));
    ASSERT_TRUE(boost::filesystem::exists(outputFiles.alphaFile));
    ASSERT_TRUE(boost::filesystem::exists(outputFiles.betaFile));
    ASSERT_TRUE(boost::filesystem::exists(outputFiles.solvationInputFile));

    std::fstream in;
    in.open(outputFiles.controlFile);
    std::string line;
    std::getline(in, line);
    bool cosmoSectionFound = line.find("$cosmo") != std::string::npos;
    ASSERT_TRUE(cosmoSectionFound);
    std::getline(in, line);
    std::getline(in, line);
    // check if correct solvent radius is applied
    bool correctSolvRadFound = line.find("rsolv= 1.93") != std::string::npos;
    ASSERT_TRUE(correctSolvRadFound);

    // Check that the states handler works.
    std::string testString = "Test";
    std::string changedContentString = "Changed Content";
    std::string content1, content2, content3;

    auto const changeContent = [](std::string& file, std::string newContent) {
      std::ofstream stream;
      stream.open(file);
      if (stream.is_open()) {
        stream << newContent << std::endl;
      }
      stream.close();
    };

    auto const checkContent = [](std::string& file, std::string& contentToCheck) {
      std::ifstream check;
      check.open(file);
      if (check.is_open()) {
        check >> contentToCheck;
      }
      check.close();
    };

    changeContent(outputFiles.alphaFile, testString);
    changeContent(outputFiles.betaFile, testString);

    auto state = calculator.getState();

    changeContent(outputFiles.alphaFile, changedContentString);
    changeContent(outputFiles.betaFile, changedContentString);

    calculator.loadState(state);

    checkContent(outputFiles.alphaFile, content1);
    checkContent(outputFiles.betaFile, content2);

    ASSERT_THAT(content1, Eq(testString));
    ASSERT_THAT(content2, Eq(testString));

    // now for the restricted case
    calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, 0);
    calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 1);
    calculator.settings().modifyString(Utils::SettingsNames::spinMode, "restricted");
    try {
      calculator.calculate("");
    }
    catch (Core::UnsuccessfulCalculationException& e) {
    }
    changeContent(outputFiles.mosFile, testString);

    state = calculator.getState();

    changeContent(outputFiles.mosFile, changedContentString);

    calculator.loadState(state);

    checkContent(outputFiles.mosFile, content3);
    ASSERT_THAT(content3, Eq(testString));

    boost::filesystem::remove_all(calculator.getCalculationDirectory());
  }
  else {
    auto logger = Core::Log();
    logger.output << "Turbomole input creation was not tested directly as no binary path was specified." << Core::Log::endl;
  }
#endif
}

TEST_F(ATurbomoleTest, PointChargeEmbeddingWorks) {
#ifndef _WIN32
  const char* envVariablePtr = std::getenv("TURBODIR");
  if (envVariablePtr) {
    auto pointChargesFile = Utils::NativeFilenames::combinePathSegments(ressourcesDir, "test_point_charges.pc");
    calculator.settings().modifyString(Scine::Utils::ExternalQC::SettingsNames::pointChargesFile, pointChargesFile);

    std::stringstream stream("3\n\n"
                             "H 3.073966 2.638248 0.173676\n"
                             "H 2.715150 1.261101 0.831928\n"
                             "O 2.828578 1.726535 0.000000\n");

    auto structure = Utils::XyzStreamHandler::read(stream);
    calculator.setStructure(structure);

    try {
      calculator.calculate("");
    }
    catch (Core::UnsuccessfulCalculationException& e) {
    }

    auto calcDir = calculator.getCalculationDirectory();
    ExternalQC::TurbomoleFiles outputFiles;
    setCorrectTurbomoleFileNames(outputFiles, calcDir);

    ASSERT_TRUE(boost::filesystem::exists(outputFiles.controlFile));

    std::fstream in;
    in.open(outputFiles.controlFile);
    std::string line;

    // check if all sections are correctly added to the control file
    bool pointChargesFound = false;
    bool pointChargesGradientFound = false;
    bool pointChargesGradientFileFound = false;

    while (std::getline(in, line)) {
      if (line.find("$point_charges") != std::string::npos)
        pointChargesFound = true;
      if (line.find("$point_charge_gradients file=pc_gradient") != std::string::npos)
        pointChargesGradientFileFound = true;
      if (line.find("$drvopt") != std::string::npos) {
        std::getline(in, line);
        if (line.find("point charges") != std::string::npos)
          pointChargesGradientFound = true;
      }
    }
    ASSERT_TRUE(pointChargesFound);
    ASSERT_TRUE(pointChargesGradientFound);
    ASSERT_TRUE(pointChargesGradientFileFound);
    boost::filesystem::remove_all(calculator.getCalculationDirectory());
  }
  else {
    auto logger = Core::Log();
    logger.output << "Turbomole input creation was not tested directly as no binary path was specified." << Core::Log::endl;
  }
#endif
}

TEST_F(ATurbomoleTest, SpinModeAssignmentIsChecked) {
#ifndef _WIN32
  const char* envVariablePtr = std::getenv("TURBODIR");
  if (envVariablePtr) {
    calculator.settings().modifyString(Utils::SettingsNames::spinMode, "restricted");

    std::stringstream stream("2\n\n"
                             "O     0.0000000    0.0000000    0.6093992\n"
                             "O     0.0000000    0.0000000    -0.6093992\n");

    auto structure = Utils::XyzStreamHandler::read(stream);
    calculator.setStructure(structure);

    ASSERT_THROW(calculator.calculate(""), std::runtime_error);

    auto calcDir = calculator.getCalculationDirectory();
    ExternalQC::TurbomoleFiles outputFiles;
    setCorrectTurbomoleFileNames(outputFiles, calcDir);

    std::fstream in;
    std::string line;
    in.open(outputFiles.controlFile);
    bool unrestrictedWasAssigned = false;
    while (std::getline(in, line)) {
      if (line.find("$uhfmo") != std::string::npos) {
        unrestrictedWasAssigned = true;
        return;
      }
    }
    ASSERT_TRUE(unrestrictedWasAssigned);

    boost::filesystem::remove_all(calculator.getCalculationDirectory());
  }
  else {
    auto logger = Core::Log();
    logger.output << "Turbomole input creation was not tested directly as no binary path was specified." << Core::Log::endl;
  }
#endif
}

TEST_F(ATurbomoleTest, BasisSetStringIsConvertedCorrectly) {
#ifndef _WIN32
  const char* envVariablePtr = std::getenv("TURBODIR");
  if (envVariablePtr) {
    auto calcDir = calculator.getCalculationDirectory();
    auto exeDir = calculator.getTurbomoleExecutableBase();
    ExternalQC::TurbomoleHelper helper(calcDir, exeDir);
    // A basis set that needs to be transformed
    calculator.settings().modifyString(Utils::SettingsNames::basisSet, "def2-svp");
    auto basisSet1 = calculator.settings().getString(Utils::SettingsNames::basisSet);
    helper.mapBasisSetToTurbomoleStringRepresentation(basisSet1);
    ASSERT_EQ(basisSet1, "def2-SVP");

    // a basis set that is upper case only
    calculator.settings().modifyString(Utils::SettingsNames::basisSet, "sto-3g");
    auto basisSet2 = calculator.settings().getString(Utils::SettingsNames::basisSet);
    helper.mapBasisSetToTurbomoleStringRepresentation(basisSet2);
    ASSERT_EQ(basisSet2, "STO-3G");

    // an unknown basis set
    calculator.settings().modifyString(Utils::SettingsNames::basisSet, "dummyBasisSet");
    auto basisSet3 = calculator.settings().getString(Utils::SettingsNames::basisSet);
    ASSERT_THROW(helper.mapBasisSetToTurbomoleStringRepresentation(basisSet3), std::runtime_error);
  }
  else {
    auto logger = Core::Log();
    logger.output << "Turbomole basis set transformation was not tested directly as no binary path was specified."
                  << Core::Log::endl;
  }
#endif
}

TEST_F(ATurbomoleTest, DFTFunctionalStringIsConvertedCorrectly) {
#ifndef _WIN32
  const char* envVariablePtr = std::getenv("TURBODIR");
  if (envVariablePtr) {
    std::stringstream stream("3\n\n"
                             "H 3.073966 2.638248 0.173676\n"
                             "H 2.715150 1.261101 0.831928\n"
                             "O 2.828578 1.726535 0.000000\n");

    auto structure = Utils::XyzStreamHandler::read(stream);
    calculator.setStructure(structure);

    auto calcDir = calculator.getCalculationDirectory();
    auto exeDir = calculator.getTurbomoleExecutableBase();
    ExternalQC::TurbomoleHelper helper(calcDir, exeDir);
    // A DFT functional that needs to be transformed
    calculator.settings().modifyString(Utils::SettingsNames::method, "b3lyp");
    auto DFTFunc = calculator.settings().getString(Utils::SettingsNames::method);
    helper.mapDftFunctionalToTurbomoleStringRepresentation(DFTFunc);
    ASSERT_EQ(DFTFunc, "b3-lyp");

    // A DFT functional that does not need to be transformed
    calculator.settings().modifyString(Utils::SettingsNames::method, "pbe");
    auto DFTFunc2 = calculator.settings().getString(Utils::SettingsNames::method);
    helper.mapDftFunctionalToTurbomoleStringRepresentation(DFTFunc2);
    ASSERT_EQ(DFTFunc2, "pbe");

    // an unknown functional should throw an error
    calculator.settings().modifyString(Utils::SettingsNames::method, "dummyDftFunc");
    auto DFTFunc3 = calculator.settings().getString(Utils::SettingsNames::method);
    ASSERT_THROW(calculator.calculate(""), std::runtime_error);

    boost::filesystem::remove_all(calculator.getCalculationDirectory());
  }
  else {
    auto logger = Core::Log();
    logger.output << "Turbomole DFT functional transformation was not tested directly as no binary path was specified."
                  << Core::Log::endl;
  }
#endif
}

TEST_F(ATurbomoleTest, ImprovedScfConvergenceSettingsAreAppliedCorrectly) {
#ifndef _WIN32
  const char* envVariablePtr = std::getenv("TURBODIR");
  if (envVariablePtr) {
    calculator.settings().modifyDouble(Utils::ExternalQC::SettingsNames::scfOrbitalShift, 0.5);
    calculator.settings().modifyBool(Utils::SettingsNames::scfDamping, true);
    calculator.settings().modifyInt(Utils::SettingsNames::maxScfIterations, 250);
    calculator.settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-4);

    auto scfIter = calculator.settings().getInt(Utils::SettingsNames::maxScfIterations);
    auto scfOrbShift = calculator.settings().getDouble(Utils::ExternalQC::SettingsNames::scfOrbitalShift);
    auto scfConvCrit = int(round(-log10(calculator.settings().getDouble(Utils::SettingsNames::selfConsistenceCriterion))));

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

    auto calcDir = calculator.getCalculationDirectory();
    ExternalQC::TurbomoleFiles outputFiles;
    setCorrectTurbomoleFileNames(outputFiles, calcDir);

    std::ifstream input;
    input.open(outputFiles.controlFile);
    auto content = std::string(std::istreambuf_iterator<char>{input}, {});
    input.close();

    std::regex scfIterRegex(R"((scfiterlimit)+ +([0-9]+))");
    std::regex scfDamp(R"(scfdamp   start=0.500  step=0.05  min=0.10)");
    std::regex scfOrbitalShift(R"((scforbitalshift closedshell=)++([-\.0-9]+))");
    std::regex convThreshold(R"((scfconv)+ +([0-9]+))");

    std::smatch m1, m2, m3, m4;
    bool scfIterFound = std::regex_search(content, m1, scfIterRegex);
    bool scfDampFound = std::regex_search(content, m2, scfDamp);
    bool scfOrbitalShiftFound = std::regex_search(content, m3, scfOrbitalShift);
    bool convThresholdFound = std::regex_search(content, m4, convThreshold);

    ASSERT_TRUE(scfIterFound);
    ASSERT_TRUE(scfDampFound);
    ASSERT_TRUE(scfOrbitalShiftFound);
    ASSERT_TRUE(convThresholdFound);

    ASSERT_EQ(std::stoi(m1[2]), scfIter);
    ASSERT_EQ(std::stod(m3[2]), scfOrbShift);
    ASSERT_EQ(std::stoi(m4[2]), scfConvCrit);
    boost::filesystem::remove_all(calculator.getCalculationDirectory());

    // Check if scf convergence criteria are adapted for Hessian and/or gradients calculation
    calculator.setRequiredProperties(Property::Energy | Property::Hessian);
    auto silentLog = Core::Log::silent();
    calculator.setLog(silentLog);
    // Trigger the applySettings() function via cloning
    auto secondCalculator = calculator.clone();
    ASSERT_THAT(secondCalculator->settings().getDouble(Utils::SettingsNames::selfConsistenceCriterion), Eq(1e-8));
    // request calculation of gradients
    calculator.setRequiredProperties(Property::Energy | Property::Gradients);
    auto thirdCalculator = calculator.clone();
    ASSERT_THAT(thirdCalculator->settings().getDouble(Utils::SettingsNames::selfConsistenceCriterion), Eq(1e-8));
  }
  else {
    auto logger = Core::Log();
    logger.output << "Turbomole control file was not checked as no binary path was specified." << Core::Log::endl;
  }
#endif
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
