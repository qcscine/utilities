/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/ExternalQC/MRCC/MrccSettings.h"
#include <Utils/ExternalQC/SettingsNames.h>
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine::Utils::ExternalQC {

MrccSettings::MrccSettings() : Settings("MrccSettings") {
  addMolecularCharge(_fields);
  addSpinMultiplicity(_fields);
  addSelfConsistenceCriterion(_fields);
  addMaxScfIterations(_fields);
  addMethod(_fields);
  addBasisSet(_fields);
  addSpinMode(_fields);
  addNumProcs(_fields);
  addMemory(_fields);
  addBaseWorkingDirectory(_fields);
  addScfDamping(_fields);
  addScfDampingValue(_fields);
  addScfOrbitalShift(_fields);
  addSolvent(_fields);
  addSolvation(_fields);
  addTemperature(_fields);
  addElectronicTemperature(_fields);
  addPressure(_fields);
  resetToDefaults();
}

void MrccSettings::addMolecularCharge(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor molecularCharge("Sets the molecular charge to use in the calculation.");
  molecularCharge.setMinimum(-10);
  molecularCharge.setMaximum(10);
  molecularCharge.setDefaultValue(0);
  settings.push_back(Scine::Utils::SettingsNames::molecularCharge, std::move(molecularCharge));
}

void MrccSettings::addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor spinMultiplicity("Sets the desired spin multiplicity to use in the "
                                                           "calculation.");
  spinMultiplicity.setMinimum(1);
  spinMultiplicity.setMaximum(10);
  spinMultiplicity.setDefaultValue(1);
  settings.push_back(Scine::Utils::SettingsNames::spinMultiplicity, std::move(spinMultiplicity));
}

void MrccSettings::addSelfConsistenceCriterion(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor selfConsistenceCriterion("Sets the desired convergence criterion.");
  selfConsistenceCriterion.setMinimum(0);
  selfConsistenceCriterion.setDefaultValue(1e-7);
  settings.push_back(Scine::Utils::SettingsNames::selfConsistenceCriterion, std::move(selfConsistenceCriterion));
}

void MrccSettings::addMethod(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor method("The method used in the MRCC calculation.");
  method.setDefaultValue("lno-ccsd(t)");
  settings.push_back(Utils::SettingsNames::method, std::move(method));
}

void MrccSettings::addBasisSet(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor basisSet("The basis set used in the calculation.");
  basisSet.setDefaultValue("def2-SVP");
  settings.push_back(Utils::SettingsNames::basisSet, std::move(basisSet));
}

void MrccSettings::addSpinMode(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor spinMode("The spin mode such as 'restricted' or 'unrestricted'.");
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Any));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Restricted));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::RestrictedOpenShell));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Unrestricted));
  spinMode.setDefaultOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Any));
  settings.push_back(Utils::SettingsNames::spinMode, std::move(spinMode));
}

void MrccSettings::addNumProcs(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor numProcs("Number of processes for the calculation.");
  numProcs.setDefaultValue(1);
  numProcs.setMinimum(1);
  settings.push_back(Utils::SettingsNames::externalProgramNProcs, std::move(numProcs));
}

void MrccSettings::addMemory(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor memory("Available system memory in MB.");
  memory.setDefaultValue(1024);
  settings.push_back(Utils::SettingsNames::externalProgramMemory, std::move(memory));
}

void MrccSettings::addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor baseWorkingDirectory("Base directory for the calculations.");
  baseWorkingDirectory.setDefaultValue(FilesystemHelpers::currentDirectory());
  settings.push_back(SettingsNames::baseWorkingDirectory, std::move(baseWorkingDirectory));
}

void MrccSettings::addMaxScfIterations(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor maxScfIterations("Maximum number of SCF iterations.");
  maxScfIterations.setMinimum(1);
  maxScfIterations.setDefaultValue(100);
  settings.push_back(Scine::Utils::SettingsNames::maxScfIterations, std::move(maxScfIterations));
}

void MrccSettings::addScfDamping(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor scfDamping("Enable stronger SCF damping (true/false).");
  scfDamping.setDefaultValue(false);
  settings.push_back(Utils::SettingsNames::scfDamping, std::move(scfDamping));
}

void MrccSettings::addScfDampingValue(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor scfDampingValue("Specify exact SCF damping value to be used.");
  scfDampingValue.setDefaultValue(0.7);
  settings.push_back(SettingsNames::scfDampingValue, std::move(scfDampingValue));
}

void MrccSettings::addScfOrbitalShift(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor scfOrbitalShift(
      "Shift the virtual orbitals to higher energies to aid convergence.");
  // internal default
  scfOrbitalShift.setDefaultValue(0.2);
  settings.push_back(SettingsNames::scfOrbitalShift, std::move(scfOrbitalShift));
}

void MrccSettings::addSolvent(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor solventOption("Sets the implicit solvent.");
  solventOption.setDefaultValue("");
  settings.push_back(Utils::SettingsNames::solvent, std::move(solventOption));
}

void MrccSettings::addSolvation(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor solvationOption(
      "Sets the implicit solvation model in the MRCC calculation.");
  solvationOption.setDefaultValue("");
  settings.push_back(Utils::SettingsNames::solvation, std::move(solvationOption));
}

void MrccSettings::addTemperature(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor value(
      "The temperature (not used by MRCC but required in the model definition).");
  value.setMinimum(0.0);
  value.setDefaultValue(298.15);
  settings.push_back(Scine::Utils::SettingsNames::temperature, std::move(value));
}

void MrccSettings::addElectronicTemperature(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor value(
      "The electronic temperature (not used by MRCC but required in the model definition).");
  value.setMinimum(0.0);
  value.setDefaultValue(0.0);
  settings.push_back(Scine::Utils::SettingsNames::electronicTemperature, std::move(value));
}

void MrccSettings::addPressure(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor value(
      "The pressure (not used by MRCC but required in the model definition).");
  value.setMinimum(0.0);
  value.setDefaultValue(101325.0);
  settings.push_back(Scine::Utils::SettingsNames::pressure, std::move(value));
}

} // namespace Scine::Utils::ExternalQC
