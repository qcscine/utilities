/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SettingPopulator.h"
#include "DescriptorCollection.h"
#include "GenericDescriptor.h"
#include <Core/Exceptions.h>

void Scine::Utils::UniversalSettings::SettingPopulator::populateLcaoSettings(SettingsCollection& settings) {
  addMolecularCharge(settings);
  addSpinMultiplicity(settings);
  addUnrestrictedCalculation(settings);
  addTemperatureOption(settings);
  addDavidsonOption(settings);
}

void Scine::Utils::UniversalSettings::SettingPopulator::populateScfSettings(SettingsCollection& settings) {
  addSelfConsistanceCriterion(settings);
  addMaxIterations(settings);
  addScfMixer(settings);
}

void Scine::Utils::UniversalSettings::SettingPopulator::populateSemiEmpiricalSettings(SettingsCollection& settings,
                                                                                      std::string defaultParameterFile) {
  Utils::UniversalSettings::FileDescriptor parameterFile("File where the parameters are stored.");
  parameterFile.setDefaultValue(std::move(defaultParameterFile));
  Utils::UniversalSettings::DirectoryDescriptor parameterRootDirectory(
      "Resource directory where all the parameters are.");
  try {
    parameterRootDirectory.setDefaultValue("");
  }
  catch (Core::InitializationException& e) {
  }

  settings.push_back(SettingsNames::parameterFile, parameterFile);
  settings.push_back(SettingsNames::parameterRootDirectory, parameterRootDirectory);
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

void Scine::Utils::UniversalSettings::SettingPopulator::addUnrestrictedCalculation(SettingsCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor unrestrictedCalculation("Switch on to run the calculation in a spin "
                                                                   "unrestricted formalism.");
  unrestrictedCalculation.setDefaultValue(false);
  settings.push_back(SettingsNames::unrestrictedCalculation, std::move(unrestrictedCalculation));
}

void Scine::Utils::UniversalSettings::SettingPopulator::addSelfConsistanceCriterion(SettingsCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor selfConsistanceCriterion("Sets the desired convergence criterion.");
  selfConsistanceCriterion.setMinimum(0.);
  selfConsistanceCriterion.setDefaultValue(1e-5);

  settings.push_back(SettingsNames::selfConsistanceCriterion, std::move(selfConsistanceCriterion));
}

void Scine::Utils::UniversalSettings::SettingPopulator::addMaxIterations(SettingsCollection& settings) {
  Utils::UniversalSettings::IntDescriptor maxIterations("Maximal number of iterations to reach self consistence.");
  maxIterations.setMinimum(1);
  maxIterations.setDefaultValue(100);

  settings.push_back(SettingsNames::maxIterations, std::move(maxIterations));
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

void Scine::Utils::UniversalSettings::SettingPopulator::addDavidsonOption(
    Scine::Utils::UniversalSettings::SettingPopulator::SettingsCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor davidsonForGroundState(
      "Use Davidson algorithm to just solve for the occupied "
      "molecular orbitals eigenpairs.");
  davidsonForGroundState.setDefaultValue(false);
  settings.push_back(SettingsNames::davidsonForGroundState, std::move(davidsonForGroundState));
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
  if (scfMixerName == SettingsNames::ScfMixers::diis)
    return Utils::scf_mixer_t::fock_diis;
  else if (scfMixerName == SettingsNames::ScfMixers::ediis)
    return Utils::scf_mixer_t::ediis;
  else if (scfMixerName == SettingsNames::ScfMixers::ediisDiis)
    return Utils::scf_mixer_t::ediis_diis;
  else if (scfMixerName == SettingsNames::ScfMixers::noMixer)
    return Utils::scf_mixer_t::none;
  else
    throw std::runtime_error("Invalid value for " + std::string(scfMixerName));
}
