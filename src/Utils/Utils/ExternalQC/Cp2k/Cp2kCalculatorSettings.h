/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_CP2KCALCULATORSETTINGS_H
#define UTILS_CP2KCALCULATORSETTINGS_H

#include <Utils/ExternalQC/SettingsNames.h>
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @class Cp2kCalculatorSettings Cp2kCalculatorSettings.h
 * @brief Settings for CP2K calculations.
 */
class Cp2kCalculatorSettings : public Scine::Utils::Settings {
 public:
  // These functions populate certain settings
  void addMolecularCharge(UniversalSettings::DescriptorCollection& settings);
  void addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings);
  void addPeriodicBoundaries(UniversalSettings::DescriptorCollection& settings);
  void addCutoff(UniversalSettings::DescriptorCollection& settings);
  void addRelCutoff(UniversalSettings::DescriptorCollection& settings);
  void addNGrids(UniversalSettings::DescriptorCollection& settings);
  void addSelfConsistenceCriterion(UniversalSettings::DescriptorCollection& settings);
  void addMaxScfIterations(UniversalSettings::DescriptorCollection& settings);
  void addMethod(UniversalSettings::DescriptorCollection& settings);
  void addBasisSet(UniversalSettings::DescriptorCollection& settings);
  void addSpinMode(UniversalSettings::DescriptorCollection& settings);
  void addNumProcs(UniversalSettings::DescriptorCollection& settings);
  void addCp2kFilenameBase(UniversalSettings::DescriptorCollection& settings);
  void addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings);
  void addDeleteTemporaryFilesOption(UniversalSettings::DescriptorCollection& settings);
  void addTemperature(UniversalSettings::DescriptorCollection& settings);
  void addPressure(UniversalSettings::DescriptorCollection& settings);
  void addScfMixing(UniversalSettings::DescriptorCollection& settings);
  void addElectronicTemperature(UniversalSettings::DescriptorCollection& settings);
  void addAdditionalMos(UniversalSettings::DescriptorCollection& settings);
  void addOrbitalTransformation(UniversalSettings::DescriptorCollection& settings);
  void addOuterScf(UniversalSettings::DescriptorCollection& settings);
  void addPoissonSolver(UniversalSettings::DescriptorCollection& settings);
  void addAllowUnconvergedScf(UniversalSettings::DescriptorCollection& settings);
  void addScfGuess(UniversalSettings::DescriptorCollection& settings);
  void addDipoleCorrection(UniversalSettings::DescriptorCollection& settings);
  void addAdditionalOutput(UniversalSettings::DescriptorCollection& settings);
  void addScfCriterionEnforce(UniversalSettings::DescriptorCollection& settings);

  /**
   * @brief Constructor that populates the Cp2kCalculatorSettings.
   */
  Cp2kCalculatorSettings() : Settings("Cp2kCalculatorSettings") {
    addMolecularCharge(_fields);
    addSpinMultiplicity(_fields);
    addSelfConsistenceCriterion(_fields);
    addPeriodicBoundaries(_fields);
    addCutoff(_fields);
    addRelCutoff(_fields);
    addNGrids(_fields);
    addMaxScfIterations(_fields);
    addMethod(_fields);
    addBasisSet(_fields);
    addSpinMode(_fields);
    addNumProcs(_fields);
    addCp2kFilenameBase(_fields);
    addBaseWorkingDirectory(_fields);
    addDeleteTemporaryFilesOption(_fields);
    addTemperature(_fields);
    addPressure(_fields);
    addScfMixing(_fields);
    addElectronicTemperature(_fields);
    addAdditionalMos(_fields);
    addOrbitalTransformation(_fields);
    addOuterScf(_fields);
    addPoissonSolver(_fields);
    addAllowUnconvergedScf(_fields);
    addScfGuess(_fields);
    addDipoleCorrection(_fields);
    addAdditionalOutput(_fields);
    addScfCriterionEnforce(_fields);
    resetToDefaults();
  };
};

inline void Cp2kCalculatorSettings::addMolecularCharge(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor molecularCharge("Sets the molecular charge to use in the calculation.");
  molecularCharge.setMinimum(-10);
  molecularCharge.setMaximum(10);
  molecularCharge.setDefaultValue(0);
  settings.push_back(Scine::Utils::SettingsNames::molecularCharge, std::move(molecularCharge));
}

inline void Cp2kCalculatorSettings::addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor spinMultiplicity("Sets the desired spin multiplicity to use in the "
                                                           "calculation.");
  spinMultiplicity.setMinimum(1);
  spinMultiplicity.setMaximum(10);
  spinMultiplicity.setDefaultValue(1);
  settings.push_back(Scine::Utils::SettingsNames::spinMultiplicity, std::move(spinMultiplicity));
}

inline void Cp2kCalculatorSettings::addCutoff(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor cutoff("Sets the plane wave cutoff of the finest grid in Ry.");
  cutoff.setMinimum(0.0);
  cutoff.setDefaultValue(300.0);
  settings.push_back(SettingsNames::planeWaveCutoff, std::move(cutoff));
}

inline void Cp2kCalculatorSettings::addRelCutoff(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor relCutoff("Determines the grid at which a Gaussian is mapped, "
                                                       "giving the cutoff in Ry used for a gaussian with alpha=1");
  relCutoff.setMinimum(0);
  relCutoff.setDefaultValue(60.0);
  settings.push_back(SettingsNames::relMultiGridCutoff, std::move(relCutoff));
}

inline void Cp2kCalculatorSettings::addNGrids(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor nGrids("Sets the desired number of grids.");
  nGrids.setMinimum(1);
  nGrids.setMaximum(10);
  nGrids.setDefaultValue(5);
  settings.push_back(SettingsNames::nGrids, std::move(nGrids));
}

inline void Cp2kCalculatorSettings::addPeriodicBoundaries(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor periodicBoundaries("Sets the unit cell.");
  periodicBoundaries.setDefaultValue(SettingsNames::moleculePeriodicBoundaries);
  settings.push_back(Utils::SettingsNames::periodicBoundaries, std::move(periodicBoundaries));
}

inline void Cp2kCalculatorSettings::addSelfConsistenceCriterion(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor selfConsistenceCriterion("Sets the desired convergence criterion.");
  selfConsistenceCriterion.setMinimum(0);
  selfConsistenceCriterion.setDefaultValue(1e-7);
  settings.push_back(Utils::SettingsNames::selfConsistenceCriterion, std::move(selfConsistenceCriterion));
}

inline void Cp2kCalculatorSettings::addMethod(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor method("The method used in the CP2K calculation.");
  method.setDefaultValue(""); // empty to avoid sneaky DFT instead of semi-empirics
  settings.push_back(Utils::SettingsNames::method, std::move(method));
}

inline void Cp2kCalculatorSettings::addBasisSet(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor basisSet("The basis set and pseudopotential used in the CP2K calculation. "
                                                      "Currently, only MOLOPT basis sets are supported.");
  basisSet.setDefaultValue("DZVP-MOLOPT-GTH");
  settings.push_back(Utils::SettingsNames::basisSet, std::move(basisSet));
}

inline void Cp2kCalculatorSettings::addSpinMode(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor spinMode("The spin mode such as 'restricted' or 'unrestricted'.");
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Any));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Restricted));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::RestrictedOpenShell));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Unrestricted));
  spinMode.setDefaultOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Any));
  settings.push_back(Utils::SettingsNames::spinMode, std::move(spinMode));
}

inline void Cp2kCalculatorSettings::addNumProcs(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor cp2kNumProcs("Number of processes for the CP2K calculation.");
  cp2kNumProcs.setDefaultValue(1);
  cp2kNumProcs.setMinimum(1);
  settings.push_back(Utils::SettingsNames::externalProgramNProcs, std::move(cp2kNumProcs));
}

inline void Cp2kCalculatorSettings::addCp2kFilenameBase(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor cp2kFilenameBase("Base of the file name of the CP2K calculations.");
  cp2kFilenameBase.setDefaultValue("cp2k_calc");
  settings.push_back(SettingsNames::cp2kFilenameBase, std::move(cp2kFilenameBase));
}

inline void Cp2kCalculatorSettings::addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor baseWorkingDirectory("Base directory for the CP2K calculations.");
  baseWorkingDirectory.setDefaultValue(FilesystemHelpers::currentDirectory());
  settings.push_back(SettingsNames::baseWorkingDirectory, std::move(baseWorkingDirectory));
}

inline void Cp2kCalculatorSettings::addMaxScfIterations(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor maxScfIterations("Maximum number of inner SCF iterations.");
  maxScfIterations.setMinimum(1);
  maxScfIterations.setDefaultValue(100);
  settings.push_back(Utils::SettingsNames::maxScfIterations, std::move(maxScfIterations));
}

inline void Cp2kCalculatorSettings::addDeleteTemporaryFilesOption(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor deleteTemporaryFiles(
      "Delete all files with the .bak extension after an CP2K calculation has failed.");
  deleteTemporaryFiles.setDefaultValue(true);
  settings.push_back(SettingsNames::deleteTemporaryFiles, std::move(deleteTemporaryFiles));
}

inline void Cp2kCalculatorSettings::addTemperature(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor temperature("Sets the temperature for the thermochemical calculation.");
  temperature.setDefaultValue(298.15);
  settings.push_back(Utils::SettingsNames::temperature, std::move(temperature));
}

inline void Cp2kCalculatorSettings::addPressure(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor pressure("Sets the pressure for the thermochemical calculation in Pa.");
  pressure.setDefaultValue(101325.0);
  settings.push_back(Utils::SettingsNames::pressure, std::move(pressure));
}

inline void Cp2kCalculatorSettings::addScfMixing(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor scfDamping("Specify SCF mixing method.");
  scfDamping.addOption("broyden_mixing");
  scfDamping.addOption("broyden_mixing_new");
  scfDamping.addOption("direct_p_mixing");
  scfDamping.addOption("kerker_mixing");
  scfDamping.addOption("multisecant_mixing");
  scfDamping.addOption("none_mixing");
  scfDamping.addOption("pulay_mixing");
  scfDamping.setDefaultOption("broyden_mixing");
  settings.push_back(Utils::SettingsNames::scfDamping, std::move(scfDamping));
}

inline void Cp2kCalculatorSettings::addElectronicTemperature(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor electronicTemperature(
      "Sets the electronic temperature for SCF calculations.");
  electronicTemperature.setMinimum(0.0);
  electronicTemperature.setDefaultValue(0.0);
  settings.push_back(Utils::SettingsNames::electronicTemperature, std::move(electronicTemperature));
}

inline void Cp2kCalculatorSettings::addAdditionalMos(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor additionalMos("Specify the number of additional molecular orbitals.");
  additionalMos.setMinimum(0);
  additionalMos.setDefaultValue(0);
  settings.push_back(SettingsNames::additionalMos, std::move(additionalMos));
}

inline void Cp2kCalculatorSettings::addOrbitalTransformation(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor ot("Specify an orbital transformation minimizer. "
                                                    "None deactivates orbital transformation.");
  ot.addOption("");
  ot.addOption("broyden");
  ot.addOption("cg");
  ot.addOption("diis");
  ot.addOption("sd");
  ot.setDefaultOption("");
  settings.push_back(SettingsNames::orbitalTransformation, std::move(ot));
}

inline void Cp2kCalculatorSettings::addOuterScf(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor outerScfIterations("Maximum number of outer SCF iterations.");
  outerScfIterations.setMinimum(0);
  outerScfIterations.setDefaultValue(0);
  settings.push_back(SettingsNames::outerScf, std::move(outerScfIterations));
}

inline void Cp2kCalculatorSettings::addPoissonSolver(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor poissonSolver(
      "Specify the poisson solver. None picks the default solver based on the periodicity.");
  poissonSolver.addOption("");
  poissonSolver.addOption("analytic");
  poissonSolver.addOption("implicit");
  poissonSolver.addOption("mt");
  poissonSolver.addOption("multipole");
  poissonSolver.addOption("periodic");
  poissonSolver.addOption("wavelet");
  poissonSolver.setDefaultOption("");
  settings.push_back(SettingsNames::poissonSolver, std::move(poissonSolver));
}

inline void Cp2kCalculatorSettings::addAllowUnconvergedScf(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor allowUnconvergedScf("Whether unconverged SCF is ignored.");
  allowUnconvergedScf.setDefaultValue(false);
  settings.push_back(SettingsNames::allowUnconvergedScf, std::move(allowUnconvergedScf));
}

inline void Cp2kCalculatorSettings::addScfGuess(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor scfGuess("The guess for the SCF. "
                                                          "Restart defaults to atomic if no restart available.");
  scfGuess.addOption("restart");
  scfGuess.addOption("atomic");
  scfGuess.addOption("core");
  scfGuess.addOption("history_restart");
  scfGuess.addOption("mopac");
  scfGuess.addOption("random");
  scfGuess.setDefaultOption("restart");
  settings.push_back(SettingsNames::scfGuess, std::move(scfGuess));
}

inline void Cp2kCalculatorSettings::addDipoleCorrection(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor dipoleCorrection("Whether a dipole correction along z-axis is applied.");
  dipoleCorrection.setDefaultValue(false);
  settings.push_back(SettingsNames::dipoleCorrection, std::move(dipoleCorrection));
}

inline void Cp2kCalculatorSettings::addAdditionalOutput(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor additionalOutput("Filename of additional output file.");
  additionalOutput.setDefaultValue("additional_output");
  settings.push_back(SettingsNames::additionalOutputFile, std::move(additionalOutput));
}

inline void Cp2kCalculatorSettings::addScfCriterionEnforce(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor scfEnforce(
      "Whether the set self_consistence_criterion should not be made stricter, "
      "even if derivative quantities are calculated.");
  scfEnforce.setDefaultValue(false);
  settings.push_back(SettingsNames::enforceScfCriterion, std::move(scfEnforce));
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_CP2KCALCULATORSETTINGS_H
