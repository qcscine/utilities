/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/ExternalQC/MRCC/MrccCCCalculator.h"
#include "Utils/ExternalQC/MRCC/MrccDFTCalculator.h"
#include "Utils/ExternalQC/MRCC/MrccHFCalculator.h"
#include "Utils/ExternalQC/MRCC/MrccHelper.h"
#include "Utils/ExternalQC/MRCC/MrccIO.h"
#include "Utils/ExternalQC/MRCC/MrccMP2Calculator.h"
#include "Utils/ExternalQC/MRCC/MrccSettings.h"
#include "Utils/IO/NativeFilenames.h"
#include <Utils/ExternalQC/SettingsNames.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <gmock/gmock.h>
#include <boost/dll/runtime_symbol_info.hpp>
#include <boost/filesystem.hpp>

using namespace testing;
namespace Scine {
namespace Utils {
namespace ExternalQC {
namespace Tests {

class MrccTest : public Test {
 public:
  std::string ressourcesDir;
  MrccSettings dftSettings;
  AtomCollection methane;

 private:
  void SetUp() final {
    dftSettings.modifyInt(Utils::SettingsNames::externalProgramNProcs, 2);
    dftSettings.modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 0.0001);
    dftSettings.modifyInt(Utils::SettingsNames::molecularCharge, 1);
    dftSettings.modifyInt(Utils::SettingsNames::spinMultiplicity, 2);
    dftSettings.modifyInt(Utils::SettingsNames::maxScfIterations, 125);
    dftSettings.modifyString(Utils::SettingsNames::method, "PBE-D3BJ");
    dftSettings.modifyString(Utils::SettingsNames::basisSet, "def2-SVP");
    dftSettings.modifyString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, "test_1");
    dftSettings.modifyString(Utils::SettingsNames::solvation, "iefpcm");
    dftSettings.modifyString(Utils::SettingsNames::solvent, "water");

    boost::filesystem::path pathToResource = boost::dll::program_location().parent_path();
    pathToResource /= "Resources";
    ressourcesDir = NativeFilenames::combinePathSegments(pathToResource.string(), "MRCC");

    std::stringstream stream2("5\n\n"
                              "C     0.00000000   0.00000001  -0.00000097\n"
                              "H     0.62612502   0.62612484   0.62613824\n"
                              "H    -0.62612503  -0.62612486   0.62613824\n"
                              "H    -0.62612481   0.62612463  -0.62613657\n"
                              "H     0.62612481  -0.62612464  -0.62613657\n");
    methane = Utils::XyzStreamHandler::read(stream2);
  }
};

TEST_F(MrccTest, SettingsAreSetCorrectly) {
  ASSERT_THAT(dftSettings.getInt(Utils::SettingsNames::externalProgramNProcs), Eq(2));
  ASSERT_THAT(dftSettings.getDouble(Utils::SettingsNames::selfConsistenceCriterion), Eq(0.0001));
  ASSERT_THAT(dftSettings.getInt(Utils::SettingsNames::molecularCharge), Eq(1));
  ASSERT_THAT(dftSettings.getInt(Utils::SettingsNames::spinMultiplicity), Eq(2));
  ASSERT_THAT(dftSettings.getInt(Utils::SettingsNames::maxScfIterations), Eq(125));
  ASSERT_THAT(dftSettings.getString(Utils::SettingsNames::method), Eq("PBE-D3BJ"));
  ASSERT_THAT(dftSettings.getString(Utils::SettingsNames::basisSet), Eq("def2-SVP"));
  ASSERT_THAT(dftSettings.getString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory), Eq("test_1"));
  ASSERT_THAT(dftSettings.getString(Utils::SettingsNames::solvation), Eq("iefpcm"));
  ASSERT_THAT(dftSettings.getString(Utils::SettingsNames::solvent), Eq("water"));
}

TEST_F(MrccTest, ReadDFTOutput) {
  MrccFiles files(ressourcesDir);
  files.output = NativeFilenames::combinePathSegments(ressourcesDir, "dft.out");
  ;
  MrccIO mrccIo(files, dftSettings, "DFT");
  const auto out = mrccIo.readOutput();
  const double energy = mrccIo.getEnergy(out);
  EXPECT_NEAR(energy, -231.782146723901, 1e-9);
}

TEST_F(MrccTest, ReadHFOutput) {
  MrccFiles files(ressourcesDir);
  auto hfSettings = dftSettings;
  hfSettings.modifyString(Utils::SettingsNames::method, "HF");
  files.output = NativeFilenames::combinePathSegments(ressourcesDir, "hf.out");
  MrccIO mrccIo(files, hfSettings, "HF");
  const auto out = mrccIo.readOutput();
  const double energy = mrccIo.getEnergy(out);
  EXPECT_NEAR(energy, -230.5359880664254888, 1e-9);
}

TEST_F(MrccTest, ReadMP2Output) {
  MrccFiles files(ressourcesDir);
  auto mp2Settings = dftSettings;
  mp2Settings.modifyString(Utils::SettingsNames::method, "LNO-MP2");
  files.output = NativeFilenames::combinePathSegments(ressourcesDir, "mp2.out");
  ;
  MrccIO mrccIo(files, mp2Settings, "MP2");
  const auto out = mrccIo.readOutput();
  const double energy = mrccIo.getEnergy(out);
  EXPECT_NEAR(energy, -231.314757631330, 1e-9);
}

TEST_F(MrccTest, ReadCCSDOutput) {
  MrccFiles files(ressourcesDir);
  auto ccsdSettings = dftSettings;
  ccsdSettings.modifyString(Utils::SettingsNames::method, "LNO-CCSD");
  files.output = NativeFilenames::combinePathSegments(ressourcesDir, "ccsd.out");
  ;
  MrccIO mrccIo(files, ccsdSettings, "CC");
  const auto out = mrccIo.readOutput();
  const double energy = mrccIo.getEnergy(out);
  EXPECT_NEAR(energy, -231.356548489232, 1e-9);
}

TEST_F(MrccTest, ReadCCSD_TOutput) {
  MrccFiles files(ressourcesDir);
  auto ccsdtSettings = dftSettings;
  ccsdtSettings.modifyString(Utils::SettingsNames::method, "LNO-CCSD(T)");
  files.output = NativeFilenames::combinePathSegments(ressourcesDir, "ccsd_t.out");
  ;
  MrccIO mrccIo(files, ccsdtSettings, "CC");
  const auto out = mrccIo.readOutput();
  const double energy = mrccIo.getEnergy(out);
  EXPECT_NEAR(energy, -231.391176960628, 1e-9);
}

TEST_F(MrccTest, runHFCalculation) {
  const char* mrccDir = std::getenv(MrccCalculator::binaryEnvVariable);
  if (mrccDir) {
    MrccHFCalculator calculator;
    calculator.settings().modifyString(Utils::SettingsNames::method, "HF");
    calculator.settings().modifyString(Utils::SettingsNames::solvation, "iefpcm");
    calculator.settings().modifyString(Utils::SettingsNames::solvent, "water");
    calculator.setRequiredProperties(Property::Energy);
    calculator.setStructure(methane);
    const auto res = calculator.calculate("");
    EXPECT_NEAR(res.get<Property::Energy>(), -40.169246, 1e-5);

    calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, 1);
    calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 2);
    const auto res2 = calculator.calculate("");
    EXPECT_NEAR(res2.get<Property::Energy>(), -39.790963, 1e-5);

    boost::filesystem::remove_all(calculator.getCalculationDirectory());
  }
}

TEST_F(MrccTest, runDFTCalculation) {
  const char* mrccDir = std::getenv(MrccCalculator::binaryEnvVariable);
  if (mrccDir) {
    MrccDFTCalculator calculator;
    calculator.settings().modifyString(Utils::SettingsNames::method, "pbe-d3bj");
    calculator.settings().modifyString(Utils::SettingsNames::solvation, "iefpcm");
    calculator.settings().modifyString(Utils::SettingsNames::solvent, "water");
    calculator.setRequiredProperties(Property::Energy);
    calculator.setStructure(methane);
    const auto res = calculator.calculate("");
    EXPECT_NEAR(res.get<Property::Energy>(), -40.415282, 1e-5);

    calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, 1);
    calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 2);
    const auto res2 = calculator.calculate("");
    EXPECT_NEAR(res2.get<Property::Energy>(), -39.99014, 1e-5);

    boost::filesystem::remove_all(calculator.getCalculationDirectory());
  }
}

TEST_F(MrccTest, runCCSDCalculation) {
  const char* mrccDir = std::getenv(MrccCalculator::binaryEnvVariable);
  if (mrccDir) {
    MrccCCCalculator calculator;
    calculator.settings().modifyString(Utils::SettingsNames::method, "LNO-CCSD");
    calculator.settings().modifyString(Utils::SettingsNames::solvation, "iefpcm");
    calculator.settings().modifyString(Utils::SettingsNames::solvent, "water");
    calculator.setRequiredProperties(Property::Energy);
    calculator.setStructure(methane);
    const auto res = calculator.calculate("");
    EXPECT_NEAR(res.get<Property::Energy>(), -40.353211, 1e-5);

    // At the moment MRCC cannot support open-shell LNO calculations.
    //    calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, 1);
    //    calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 2);
    //    const auto res2 = calculator.calculate("");

    boost::filesystem::remove_all(calculator.getCalculationDirectory());
  }
}

TEST_F(MrccTest, runCCSD_TCalculation) {
  const char* mrccDir = std::getenv(MrccCalculator::binaryEnvVariable);
  if (mrccDir) {
    MrccCCCalculator calculator;
    calculator.settings().modifyString(Utils::SettingsNames::method, "LNO-CCSD(T)");
    calculator.settings().modifyString(Utils::SettingsNames::solvation, "iefpcm");
    calculator.settings().modifyString(Utils::SettingsNames::solvent, "water");
    calculator.setRequiredProperties(Property::Energy);
    calculator.setStructure(methane);
    const auto res = calculator.calculate("");
    EXPECT_NEAR(res.get<Property::Energy>(), -40.356931, 1e-5);

    // At the moment MRCC cannot support open-shell LNO calculations.
    //    calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, 1);
    //    calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 2);
    //    const auto res2 = calculator.calculate("");

    boost::filesystem::remove_all(calculator.getCalculationDirectory());
  }
}

TEST_F(MrccTest, runMP2Calculation) {
  const char* mrccDir = std::getenv(MrccCalculator::binaryEnvVariable);
  if (mrccDir) {
    MrccMP2Calculator calculator;
    calculator.settings().modifyString(Utils::SettingsNames::method, "LNO-MP2");
    calculator.settings().modifyString(Utils::SettingsNames::solvation, "iefpcm");
    calculator.settings().modifyString(Utils::SettingsNames::solvent, "water");
    calculator.setRequiredProperties(Property::Energy);
    calculator.setStructure(methane);
    const auto res = calculator.calculate("");
    EXPECT_NEAR(res.get<Property::Energy>(), -40.330257, 1e-5);

    // At the moment MRCC cannot support open-shell LNO calculations.
    //    calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, 1);
    //    calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 2);
    //    const auto res2 = calculator.calculate("");

    boost::filesystem::remove_all(calculator.getCalculationDirectory());
  }
}

} // namespace Tests
} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
