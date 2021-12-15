/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SettingPopulator.h"
#include "DescriptorCollection.h"
#include "GenericDescriptor.h"
#include "Utils/Scf/LcaoUtils/SpinMode.h"
#include <Core/Exceptions.h>

void Scine::Utils::UniversalSettings::SettingPopulator::addLogOption(
    Scine::Utils::UniversalSettings::SettingPopulator::SettingsCollection& settings) {
  Utils::UniversalSettings::StringDescriptor logVerbosity("Sets the verbosity of the logger.");
  logVerbosity.setDefaultValue("output");
  settings.push_back(Utils::SettingsNames::loggerVerbosity, std::move(logVerbosity));
}

void Scine::Utils::UniversalSettings::SettingPopulator::populateLcaoSettings(SettingsCollection& settings) {
  addMolecularCharge(settings);
  addSpinMultiplicity(settings);
  addSpinMode(settings);
  addTemperatureOption(settings);
  addElectronicTemperatureOption(settings);
  addSymmetryNumberOption(settings);
}

void Scine::Utils::UniversalSettings::SettingPopulator::populateScfSettings(SettingsCollection& settings) {
  addSelfConsistenceCriterion(settings);
  addMaxScfIterations(settings);
  addScfMixer(settings);
}

void Scine::Utils::UniversalSettings::SettingPopulator::populateSemiEmpiricalSettings(SettingsCollection& settings,
                                                                                      std::string defaultParameterFile) {
  Utils::UniversalSettings::FileDescriptor parameterPath("Filesystem path where method parameters are stored.");
  parameterPath.setDefaultValue(std::move(defaultParameterFile));
  settings.push_back(SettingsNames::methodParameters, parameterPath);
}

void Scine::Utils::UniversalSettings::SettingPopulator::addMolecularCharge(SettingsCollection& settings) {
  Utils::UniversalSettings::IntDescriptor molecularCharge("Sets the molecular charge to use in the calculation.");
  molecularCharge.setMinimum(-20);
  molecularCharge.setMaximum(+20);
  molecularCharge.setDefaultValue(0);

  settings.push_back(SettingsNames::molecularCharge, std::move(molecularCharge));
}

void Scine::Utils::UniversalSettings::SettingPopulator::addSpinMultiplicity(SettingsCollection& settings) {
  Utils::UniversalSettings::IntDescriptor spinMultiplicity("Sets the desired spin multiplicity to use in the "
                                                           "calculation.");
  spinMultiplicity.setMinimum(1);
  spinMultiplicity.setMaximum(10);
  spinMultiplicity.setDefaultValue(1);

  settings.push_back(SettingsNames::spinMultiplicity, std::move(spinMultiplicity));
}

void Scine::Utils::UniversalSettings::SettingPopulator::addSpinMode(SettingsCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor spinMode("Run the calculation in a restricted or "
                                                          "unrestricted spin formalism.");
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Any));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Restricted));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::RestrictedOpenShell));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Unrestricted));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::None));
  spinMode.setDefaultOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Any));
  settings.push_back(SettingsNames::spinMode, std::move(spinMode));
}

void Scine::Utils::UniversalSettings::SettingPopulator::addSelfConsistenceCriterion(SettingsCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor scfEnergyThreshold(
      "Sets the desired convergence criterion for the energy difference.");
  scfEnergyThreshold.setMinimum(0.);
  scfEnergyThreshold.setDefaultValue(1e-7);

  settings.push_back(SettingsNames::selfConsistenceCriterion, std::move(scfEnergyThreshold));

  Utils::UniversalSettings::DoubleDescriptor densityRmsd(
      "Sets the desired convergence criterion for the density matrix RMSD.");
  densityRmsd.setMinimum(0.);
  densityRmsd.setDefaultValue(1e-5);

  settings.push_back(SettingsNames::densityRmsdCriterion, std::move(densityRmsd));
}

void Scine::Utils::UniversalSettings::SettingPopulator::addMaxScfIterations(SettingsCollection& settings) {
  Utils::UniversalSettings::IntDescriptor maxScfIterations("Maximal number of iterations to reach self consistence.");
  maxScfIterations.setMinimum(1);
  maxScfIterations.setDefaultValue(100);

  settings.push_back(SettingsNames::maxScfIterations, std::move(maxScfIterations));
}

void Scine::Utils::UniversalSettings::SettingPopulator::addScfMixer(SettingsCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor mixer("Convergence acceleration to use.");
  mixer.addOption(SettingsNames::ScfMixers::noMixer);
  mixer.addOption(SettingsNames::ScfMixers::diis);
  mixer.addOption(SettingsNames::ScfMixers::ediis);
  mixer.addOption(SettingsNames::ScfMixers::ediisDiis);
  mixer.setDefaultOption(SettingsNames::ScfMixers::diis);

  settings.push_back(SettingsNames::mixer, mixer);
}

void Scine::Utils::UniversalSettings::SettingPopulator::addTemperatureOption(
    Scine::Utils::UniversalSettings::SettingPopulator::SettingsCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor temperature("Temperature to use for thermochemical calculation.");
  temperature.setMinimum(0.);
  temperature.setMaximum(10000);
  temperature.setDefaultValue(298.15);

  settings.push_back(SettingsNames::temperature, std::move(temperature));
}

void Scine::Utils::UniversalSettings::SettingPopulator::addElectronicTemperatureOption(
    Scine::Utils::UniversalSettings::SettingPopulator::SettingsCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor etemperature("Temperature to use in the electron distribution.");
  etemperature.setMinimum(0.);
  etemperature.setMaximum(0);
  etemperature.setDefaultValue(0.0);

  settings.push_back(SettingsNames::electronicTemperature, std::move(etemperature));
}

void Scine::Utils::UniversalSettings::SettingPopulator::addSymmetryNumberOption(
    Scine::Utils::UniversalSettings::SettingPopulator::SettingsCollection& settings) {
  Utils::UniversalSettings::IntDescriptor symmetryNumber(
      "Molecular symmetry number to use for thermochemical calculation.");
  symmetryNumber.setMinimum(1);
  symmetryNumber.setDefaultValue(1);
  settings.push_back(SettingsNames::symmetryNumber, std::move(symmetryNumber));
}

std::string Scine::Utils::UniversalSettings::SettingPopulator::scfMixerToString(Utils::scf_mixer_t mixer) {
  switch (mixer) {
    case Utils::scf_mixer_t::none:
      return SettingsNames::ScfMixers::noMixer;
    case Utils::scf_mixer_t::fock_diis:
      return SettingsNames::ScfMixers::diis;
    case Utils::scf_mixer_t::ediis:
      return SettingsNames::ScfMixers::ediis;
    case Utils::scf_mixer_t::ediis_diis:
      return SettingsNames::ScfMixers::ediisDiis;
    default:
      throw std::runtime_error("Unknown conversion from Utils::scf_mixer_t to std::string. Enum id is " +
                               std::to_string(static_cast<int>(mixer)));
  }
}

Scine::Utils::scf_mixer_t Scine::Utils::UniversalSettings::SettingPopulator::stringToScfMixer(const std::string& scfMixerName) {
  if (scfMixerName == SettingsNames::ScfMixers::diis) {
    return Utils::scf_mixer_t::fock_diis;
  }
  if (scfMixerName == SettingsNames::ScfMixers::ediis) {
    return Utils::scf_mixer_t::ediis;
  }
  if (scfMixerName == SettingsNames::ScfMixers::ediisDiis) {
    return Utils::scf_mixer_t::ediis_diis;
  }
  if (scfMixerName == SettingsNames::ScfMixers::noMixer) {
    return Utils::scf_mixer_t::none;
  }

  throw std::runtime_error("Invalid value for " + std::string(scfMixerName));
}
