/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_GAUSSIANCALCULATORSETTINGS_H
#define UTILS_GAUSSIANCALCULATORSETTINGS_H

#include <Utils/ExternalQC/SettingsNames.h>
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @class GaussianCalculatorSettings GaussianCalculatorSettings.h
 * @brief Settings for Gaussian calculations.
 */
class GaussianCalculatorSettings : public Scine::Utils::Settings {
 public:
  // These functions populate certain settings
  void addMolecularCharge(UniversalSettings::DescriptorCollection& settings);
  void addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings);
  void addSelfConsistenceCriterion(UniversalSettings::DescriptorCollection& settings);
  void addMethod(UniversalSettings::DescriptorCollection& settings);
  void addGaussianFilenameBase(UniversalSettings::DescriptorCollection& settings);
  void addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings);
  void addNumProcs(UniversalSettings::DescriptorCollection& settings);
  void addBasisSet(UniversalSettings::DescriptorCollection& settings);
  void addSpinMode(UniversalSettings::DescriptorCollection& settings);
  void addMemoryForGaussian(UniversalSettings::DescriptorCollection& settings);
  void addSolvent(UniversalSettings::DescriptorCollection& settings);
  void addSolvation(UniversalSettings::DescriptorCollection& settings);
  void addElectronicTemperature(UniversalSettings::DescriptorCollection& settings);
  void addScfGuess(UniversalSettings::DescriptorCollection& settings);

  /**
   * @brief Constructor that populates the GaussianCalculatorSettings.
   */
  GaussianCalculatorSettings() : Settings("GaussianCalculatorSettings") {
    addMolecularCharge(_fields);
    addSpinMultiplicity(_fields);
    addSelfConsistenceCriterion(_fields);
    addMethod(_fields);
    addBasisSet(_fields);
    addSpinMode(_fields);
    addGaussianFilenameBase(_fields);
    addBaseWorkingDirectory(_fields);
    addNumProcs(_fields);
    addMemoryForGaussian(_fields);
    addSolvent(_fields);
    addSolvation(_fields);
    addElectronicTemperature(_fields);
    addScfGuess(_fields);
    resetToDefaults();
  };
};

inline void GaussianCalculatorSettings::addMolecularCharge(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor molecularCharge("Sets the molecular charge to use in the calculation.");
  molecularCharge.setMinimum(-10);
  molecularCharge.setMaximum(10);
  molecularCharge.setDefaultValue(0);
  settings.push_back(Scine::Utils::SettingsNames::molecularCharge, std::move(molecularCharge));
}

inline void GaussianCalculatorSettings::addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor spinMultiplicity("Sets the desired spin multiplicity to use in the "
                                                           "calculation.");
  spinMultiplicity.setMinimum(1);
  spinMultiplicity.setMaximum(10);
  spinMultiplicity.setDefaultValue(1);
  settings.push_back(Scine::Utils::SettingsNames::spinMultiplicity, std::move(spinMultiplicity));
}

inline void GaussianCalculatorSettings::addSelfConsistenceCriterion(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor selfConsistenceCriterion("Sets the desired convergence criterion.");
  selfConsistenceCriterion.setMinimum(0);
  selfConsistenceCriterion.setMaximum(1);
  selfConsistenceCriterion.setDefaultValue(1e-7);
  settings.push_back(Scine::Utils::SettingsNames::selfConsistenceCriterion, std::move(selfConsistenceCriterion));
}

inline void GaussianCalculatorSettings::addMethod(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor method("The method used in the Gaussian calculation.");
  method.setDefaultValue("PBEPBE");
  settings.push_back(Utils::SettingsNames::method, std::move(method));
}

inline void GaussianCalculatorSettings::addBasisSet(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor basisSet("The basis set used in the Gaussian calculation.");
  basisSet.setDefaultValue("def2SVP");
  settings.push_back(Utils::SettingsNames::basisSet, std::move(basisSet));
}

inline void GaussianCalculatorSettings::addSpinMode(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor spinMode("The spin mode such as 'restricted' or 'unrestricted'.");
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Any));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Restricted));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::RestrictedOpenShell));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Unrestricted));
  spinMode.setDefaultOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Any));
  settings.push_back(Utils::SettingsNames::spinMode, std::move(spinMode));
}

inline void GaussianCalculatorSettings::addGaussianFilenameBase(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor gaussianFilenameBase(
      "Base of the file name of the Gaussian calculations.");
  gaussianFilenameBase.setDefaultValue("gaussian_calc");
  settings.push_back(SettingsNames::gaussianFilenameBase, std::move(gaussianFilenameBase));
}

inline void GaussianCalculatorSettings::addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor baseWorkingDirectory("Base directory for the Gaussian calculations.");
  baseWorkingDirectory.setDefaultValue(FilesystemHelpers::currentDirectory());
  settings.push_back(SettingsNames::baseWorkingDirectory, std::move(baseWorkingDirectory));
}

inline void GaussianCalculatorSettings::addNumProcs(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor gaussianNumProcs("Number of processes for the Gaussian calculation.");
  gaussianNumProcs.setDefaultValue(1);
  gaussianNumProcs.setMinimum(1);
  settings.push_back(Utils::SettingsNames::externalProgramNProcs, std::move(gaussianNumProcs));
}

inline void GaussianCalculatorSettings::addMemoryForGaussian(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor gaussianMemory("Memory that can be used by the Gaussian calculation.");
  gaussianMemory.setDefaultValue(1024);
  settings.push_back(Utils::SettingsNames::externalProgramMemory, std::move(gaussianMemory));
}

inline void GaussianCalculatorSettings::addSolvent(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor solventOption(
      "Sets the implicit solvent using the CPCM model to be applied in the Gaussian calculation.");
  solventOption.setDefaultValue("");
  settings.push_back(Utils::SettingsNames::solvent, std::move(solventOption));
}

inline void GaussianCalculatorSettings::addSolvation(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor solvationOption(
      "Sets the implicit solvent model in the Gaussian calculation.");
  solvationOption.setDefaultValue("");
  settings.push_back(Utils::SettingsNames::solvation, std::move(solvationOption));
}

inline void GaussianCalculatorSettings::addElectronicTemperature(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor electronicTemperature(
      "Sets the electronic temperature for SCF calculations.");
  electronicTemperature.setMinimum(0.0);
  electronicTemperature.setDefaultValue(0.0);
  settings.push_back(Utils::SettingsNames::electronicTemperature, std::move(electronicTemperature));
}

inline void GaussianCalculatorSettings::addScfGuess(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor scfGuess(
      "The guess for the SCF."
      "See https://gaussian.com/guess/ for further information."
      "Read defaults to harris if no previous solution is present.");
  // There are more options in Gaussian that could be added here. See: https://gaussian.com/guess/
  scfGuess.addOption("read"); // Reads in chk file from previous run if present
  scfGuess.addOption("harris");
  scfGuess.addOption("huckel");
  scfGuess.addOption("core");
  scfGuess.addOption("only");
  scfGuess.addOption("(only, read)"); // Reads in old orbitals and calculates properties from these
  scfGuess.setDefaultOption("read");
  settings.push_back(SettingsNames::scfGuess, std::move(scfGuess));
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_GAUSSIANCALCULATORSETTINGS_H
