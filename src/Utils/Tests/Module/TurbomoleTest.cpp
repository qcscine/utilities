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
#include <Utils/ExternalQC/Turbomole/TurbomoleInputFileCreator.h>
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
  calculator.settings().modifyString(Utils::SettingsNames::method, "PBE D3BJ");
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
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::method), Eq("PBE D3BJ"));
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
    ASSERT_THROW(calculator.calculate(""), std::runtime_error);
    calculator.settings().modifyString(Utils::SettingsNames::method, "pbe d3bj");

    // charge and multiplicity pair does not match
    calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 2);
    calculator.settings().modifyString(Utils::SettingsNames::spinMode, "unrestricted");
    calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, 0);
    ASSERT_THROW(calculator.calculate(""), std::runtime_error);
    calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 1);

    // spin mode is given incorrectly
    calculator.settings().modifyString(Utils::SettingsNames::spinMode, "restricted_open_shell");
    ASSERT_THROW(calculator.calculate(""), std::runtime_error);
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
    ASSERT_THAT(results.get<Property::Energy>(), DoubleNear(-40.41523068099, 1e-8));

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
}

TEST_F(ATurbomoleTest, InputFileIsWrittenCorrectly) {
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

TEST_F(ATurbomoleTest, improvedScfConvergenceSettingsAreAppliedCorrectly) {
#ifndef _WIN32
  const char* envVariablePtr = std::getenv("TURBODIR");
  if (envVariablePtr) {
    calculator.settings().modifyDouble(Utils::ExternalQC::SettingsNames::scfOrbitalShift, 0.5);
    calculator.settings().modifyBool(Utils::SettingsNames::scfDamping, true);
    calculator.settings().modifyInt(Utils::SettingsNames::maxScfIterations, 250);

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
    std::regex scfDamp(R"(scfdamp   start=8.500  step=0.10  min=0.10)");
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
