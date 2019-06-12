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

void Scine::Utils::UniversalSettings::SettingPopulator::populateLCAOSettings(SettingsCollection& settings) {
  addMolecularCharge(settings);
  addSpinMultiplicity(settings);
  addUnrestrictedCalculation(settings);
  addLogOption(settings);
}

void Scine::Utils::UniversalSettings::SettingPopulator::populateSCFSettings(SettingsCollection& settings) {
  addSelfConsistanceCriterion(settings);
  addMaxIterations(settings);
  addSCFMixer(settings);
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
  molecularCharge.setMinimum(-10);
  molecularCharge.setMaximum(+10);
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

void Scine::Utils::UniversalSettings::SettingPopulator::addLogOption(SettingsCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor loggerVerbosity("Level of verbosity.");
  loggerVerbosity.addOption(SettingsNames::LogLevels::none);
  loggerVerbosity.addOption(SettingsNames::LogLevels::trace);
  loggerVerbosity.addOption(SettingsNames::LogLevels::debug);
  loggerVerbosity.addOption(SettingsNames::LogLevels::info);
  loggerVerbosity.addOption(SettingsNames::LogLevels::warning);
  loggerVerbosity.addOption(SettingsNames::LogLevels::error);
  loggerVerbosity.addOption(SettingsNames::LogLevels::fatal);
  loggerVerbosity.setDefaultOption(SettingsNames::LogLevels::info);

  settings.push_back(SettingsNames::loggerVerbosity, std::move(loggerVerbosity));
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

void Scine::Utils::UniversalSettings::SettingPopulator::addSCFMixer(SettingsCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor mixer("Convergence acceleration to use.");
  mixer.addOption(SettingsNames::SCFMixers::noMixer);
  mixer.addOption(SettingsNames::SCFMixers::diis);
  mixer.addOption(SettingsNames::SCFMixers::ediis);
  mixer.addOption(SettingsNames::SCFMixers::ediisDiis);
  mixer.setDefaultOption(SettingsNames::SCFMixers::diis);

  settings.push_back(SettingsNames::mixer, mixer);
}

std::string Scine::Utils::UniversalSettings::SettingPopulator::scfMixerToString(Utils::scf_mixer_t mixer) {
  switch (mixer) {
    case Utils::scf_mixer_t::none:
      return SettingsNames::SCFMixers::noMixer;
    case Utils::scf_mixer_t::fock_diis:
      return SettingsNames::SCFMixers::diis;
    case Utils::scf_mixer_t::ediis:
      return SettingsNames::SCFMixers::ediis;
    case Utils::scf_mixer_t::ediis_diis:
      return SettingsNames::SCFMixers::ediisDiis;
    default:
      throw std::runtime_error("Unknown conversion from Utils::scf_mixer_t to std::string. Enum id is " +
                               std::to_string(static_cast<int>(mixer)));
  }
}

Scine::Utils::scf_mixer_t Scine::Utils::UniversalSettings::SettingPopulator::stringToSCFMixer(const std::string& scfMixerName) {
  if (scfMixerName == SettingsNames::SCFMixers::diis)
    return Utils::scf_mixer_t::fock_diis;
  else if (scfMixerName == SettingsNames::SCFMixers::ediis)
    return Utils::scf_mixer_t::ediis;
  else if (scfMixerName == SettingsNames::SCFMixers::ediisDiis)
    return Utils::scf_mixer_t::ediis_diis;
  else if (scfMixerName == SettingsNames::SCFMixers::noMixer)
    return Utils::scf_mixer_t::none;
  else
    throw std::runtime_error("Invalid value for " + std::string(scfMixerName));
}

Scine::Utils::Log::severity_level
Scine::Utils::UniversalSettings::SettingPopulator::stringToLogLevel(const std::string& logVerbosityString) {
  if (logVerbosityString == Utils::SettingsNames::LogLevels::trace)
    return Utils::Log::severity_level::trace;
  else if (logVerbosityString == Utils::SettingsNames::LogLevels::debug)
    return Utils::Log::severity_level::debug;
  else if (logVerbosityString == Utils::SettingsNames::LogLevels::info)
    return Utils::Log::severity_level::info;
  else if (logVerbosityString == Utils::SettingsNames::LogLevels::warning)
    return Utils::Log::severity_level::warning;
  else if (logVerbosityString == Utils::SettingsNames::LogLevels::error)
    return Utils::Log::severity_level::error;
  else if (logVerbosityString == Utils::SettingsNames::LogLevels::fatal)
    return Utils::Log::severity_level::fatal;
  else
    throw Utils::UniversalSettings::OptionDoesNotExistException(logVerbosityString, Utils::SettingsNames::loggerVerbosity);
}
