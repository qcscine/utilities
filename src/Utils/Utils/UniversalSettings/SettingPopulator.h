/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SETTINGPOPULATOR_H
#define UTILS_SETTINGPOPULATOR_H

#include "SettingsNames.h"
#include <Utils/IO/Logger.h>
#include <Utils/MethodEssentials/MethodFactories/MixerFactory.h>
#include <string>
namespace Scine {
namespace Utils {
namespace UniversalSettings {

class DescriptorCollection;

/**
 * @class SettingPopulator SettingPopulator.h
 * @brief This class populates the common settings of many calculators.
 * These settings include molecular charge, spin multiplicity, restricted/unrestricted formalism, and SCF options.
 * It populates a Utils::UniversalSettings::DescriptorCollection with default values.
 */
class SettingPopulator {
 public:
  using SettingsCollection = Utils::UniversalSettings::DescriptorCollection;

  //! Populate Settings in the by-reference argument with default settings common to all LCAO methods.
  static void populateLCAOSettings(SettingsCollection& settings);
  //! Populate Settings in the by-reference argument with default settings common to all SCF methods.
  static void populateSCFSettings(SettingsCollection& settings);
  /** @brief Populate Settings in the by-reference argument with default settings common to all SemiEmpirical methods.
   *  Here the parameter files settings are set. In order for this to work SPARROWInitializer::initialize()
   *  must already be called.
   *  @param defaultParameterFile describes the default name of the parameter for a semi-empirical method.
   *         For NDDO methods it is usually parameter.xml, for DFTB it is different.
   */
  static void populateSemiEmpiricalSettings(SettingsCollection& settings,
                                            std::string defaultParameterFile = "parameter.xml");
  /**
   * @brief This function converts a mixer name into a mixer type.
   * @param scfMixerName Name of the mixer type.
   * @return a Util::scf_mixer_t, the type of a scf mixer.
   */
  static Utils::scf_mixer_t stringToSCFMixer(const std::string& scfMixerName);
  /**
   * @brief This function converts a mixer type into a mixer name.
   * @param scfMixerType Type of the scf mixer.
   * @return a std::string containing the name of the scf mixer.
   */
  static std::string scfMixerToString(Utils::scf_mixer_t mixerType);
  /**
   * @brief This function converts a log verbosity string into a severity_level type.
   * @param logVerbosityString Name of the log verbosity level.
   * @return Utils::Log::severity_level The severity_level type corresponding to the given string.
   */
  static Utils::Log::severity_level stringToLogLevel(const std::string& logVerbosityString);

 private:
  // LCAO Settings
  static void addMolecularCharge(SettingsCollection& settings);
  static void addSpinMultiplicity(SettingsCollection& settings);
  static void addUnrestrictedCalculation(SettingsCollection& settings);
  static void addLogOption(SettingsCollection& settings);

  // SCF Settings
  static void addSelfConsistanceCriterion(SettingsCollection& settings);
  static void addMaxIterations(SettingsCollection& settings);
  static void addSCFMixer(SettingsCollection& settings);
};

} // namespace UniversalSettings
} // namespace Utils
} // namespace Scine

#endif // UTILS_SETTINGPOPULATOR_H
