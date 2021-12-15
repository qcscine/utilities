/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_TURBOMOLECALCULATORSETTINGS_H
#define UTILS_TURBOMOLECALCULATORSETTINGS_H

#include <Utils/ExternalQC/SettingsNames.h>
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @class TurbomoleCalculatorSettings TurbomoleCalculatorSettings.h
 * @brief Settings for Turbomole calculations.
 */
class TurbomoleCalculatorSettings : public Scine::Utils::Settings {
 public:
  // These functions populate certain settings
  void addMolecularCharge(UniversalSettings::DescriptorCollection& settings);
  void addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings);
  void addSelfConsistenceCriterion(UniversalSettings::DescriptorCollection& settings);
  void addMaxScfIterations(UniversalSettings::DescriptorCollection& settings);
  void addMethod(UniversalSettings::DescriptorCollection& settings);
  void addBasisSet(UniversalSettings::DescriptorCollection& settings);
  void addSpinMode(UniversalSettings::DescriptorCollection& settings);
  void addNumProcs(UniversalSettings::DescriptorCollection& settings);
  void addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings);
  void addTemperature(UniversalSettings::DescriptorCollection& settings);
  void addElectronicTemperature(UniversalSettings::DescriptorCollection& settings);
  void addScfDamping(UniversalSettings::DescriptorCollection& settings);
  void addScfOrbitalShift(UniversalSettings::DescriptorCollection& settings);
  void addSolvent(UniversalSettings::DescriptorCollection& settings);
  void addSolvation(UniversalSettings::DescriptorCollection& settings);
  void addSteerOrbitals(UniversalSettings::DescriptorCollection& settings);

  /**
   * @brief Constructor that populates the TurbomoleCalculatorSettings.
   */
  TurbomoleCalculatorSettings() : Settings("TurbomoleCalculatorSettings") {
    addMolecularCharge(_fields);
    addSpinMultiplicity(_fields);
    addSelfConsistenceCriterion(_fields);
    addMaxScfIterations(_fields);
    addMethod(_fields);
    addBasisSet(_fields);
    addSpinMode(_fields);
    addNumProcs(_fields);
    addBaseWorkingDirectory(_fields);
    addTemperature(_fields);
    addScfDamping(_fields);
    addScfOrbitalShift(_fields);
    addElectronicTemperature(_fields);
    addSolvent(_fields);
    addSolvation(_fields);
    addSteerOrbitals(_fields);
    resetToDefaults();
  };
};

inline void TurbomoleCalculatorSettings::addMolecularCharge(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor molecularCharge("Sets the molecular charge to use in the calculation.");
  molecularCharge.setMinimum(-10);
  molecularCharge.setMaximum(10);
  molecularCharge.setDefaultValue(0);
  settings.push_back(Scine::Utils::SettingsNames::molecularCharge, std::move(molecularCharge));
}

inline void TurbomoleCalculatorSettings::addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor spinMultiplicity("Sets the desired spin multiplicity to use in the "
                                                           "calculation.");
  spinMultiplicity.setMinimum(1);
  spinMultiplicity.setMaximum(10);
  spinMultiplicity.setDefaultValue(1);
  settings.push_back(Scine::Utils::SettingsNames::spinMultiplicity, std::move(spinMultiplicity));
}

inline void TurbomoleCalculatorSettings::addSelfConsistenceCriterion(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor selfConsistenceCriterion("Sets the desired convergence criterion.");
  selfConsistenceCriterion.setMinimum(0);
  selfConsistenceCriterion.setDefaultValue(1e-7);
  settings.push_back(Scine::Utils::SettingsNames::selfConsistenceCriterion, std::move(selfConsistenceCriterion));
}

inline void TurbomoleCalculatorSettings::addMethod(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor method("The method used in the Turbomole calculation.");
  method.setDefaultValue("pbe");
  settings.push_back(Utils::SettingsNames::method, std::move(method));
}

inline void TurbomoleCalculatorSettings::addBasisSet(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor basisSet("The basis set used in the Turbomole calculation.");
  basisSet.setDefaultValue("def-SV(P)"); // Turbomole's internal default
  settings.push_back(Utils::SettingsNames::basisSet, std::move(basisSet));
}

inline void TurbomoleCalculatorSettings::addSpinMode(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor spinMode("The spin mode such as 'restricted' or 'unrestricted'.");
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Any));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Restricted));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::RestrictedOpenShell));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Unrestricted));
  spinMode.setDefaultOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Any));
  settings.push_back(Utils::SettingsNames::spinMode, std::move(spinMode));
}

inline void TurbomoleCalculatorSettings::addNumProcs(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor turbomoleNumProcs("Number of processes for the Turbomole calculation.");
  turbomoleNumProcs.setDefaultValue(1);
  turbomoleNumProcs.setMinimum(1);
  settings.push_back(Utils::SettingsNames::externalProgramNProcs, std::move(turbomoleNumProcs));
}

inline void TurbomoleCalculatorSettings::addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor baseWorkingDirectory("Base directory for the Turbomole calculations.");
  baseWorkingDirectory.setDefaultValue(FilesystemHelpers::currentDirectory());
  settings.push_back(SettingsNames::baseWorkingDirectory, std::move(baseWorkingDirectory));
}

inline void TurbomoleCalculatorSettings::addMaxScfIterations(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor maxScfIterations("Maximum number of SCF iterations.");
  maxScfIterations.setMinimum(1);
  maxScfIterations.setDefaultValue(100);
  settings.push_back(Scine::Utils::SettingsNames::maxScfIterations, std::move(maxScfIterations));
}

inline void TurbomoleCalculatorSettings::addTemperature(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor temperature("Sets the temperature for the thermochemical calculation.");
  temperature.setDefaultValue(298.15);
  settings.push_back(Utils::SettingsNames::temperature, std::move(temperature));
}

inline void TurbomoleCalculatorSettings::addScfDamping(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor scfDamping("Switch SCF damping on/off.");
  scfDamping.setDefaultValue(false);
  settings.push_back(Utils::SettingsNames::scfDamping, std::move(scfDamping));
}

inline void TurbomoleCalculatorSettings::addScfOrbitalShift(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor scfOrbitalShift(
      "Shift closed shells to lower energies to aid convergence.");
  // internal default
  scfOrbitalShift.setDefaultValue(0.1);
  settings.push_back(SettingsNames::scfOrbitalShift, std::move(scfOrbitalShift));
}

inline void TurbomoleCalculatorSettings::addElectronicTemperature(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor electronicTemperature(
      "Sets the electronic temperature for SCF calculations.");
  electronicTemperature.setMinimum(0.0);
  electronicTemperature.setDefaultValue(0.0);
  settings.push_back(Utils::SettingsNames::electronicTemperature, std::move(electronicTemperature));
}

inline void TurbomoleCalculatorSettings::addSolvent(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor solventOption("Sets the implicit solvent.");
  solventOption.setDefaultValue("");
  settings.push_back(Utils::SettingsNames::solvent, std::move(solventOption));
}

inline void TurbomoleCalculatorSettings::addSolvation(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor solvationOption(
      "Sets the implicit solvation model in the TURBOMOLE calculation. Currently,nothing is available.");
  solvationOption.setDefaultValue("");
  settings.push_back(Utils::SettingsNames::solvation, std::move(solvationOption));
}

inline void TurbomoleCalculatorSettings::addSteerOrbitals(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor steerOrbitals(
      "Converts internal coordinates used by default to cartesian coordinates.");
  steerOrbitals.setDefaultValue(false);
  settings.push_back(SettingsNames::steerOrbitals, std::move(steerOrbitals));
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_TURBOMOLECALCULATORSETTINGS_H
