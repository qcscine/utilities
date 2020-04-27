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
#include <Utils/Scf/ConvergenceAccelerators/ConvergenceAcceleratorFactory.h>
#include <string>
namespace Scine {
namespace Utils {
namespace UniversalSettings {

class DescriptorCollection;

/**
 * @class SettingPopulator SettingPopulator.h
 * @brief This class populates the common settings of many calculators.
 * These settings include molecular charge, spin multiplicity, restricted/unrestricted formalism, and Scf options.
 * It populates a Utils::UniversalSettings::DescriptorCollection with default values.
 */
class SettingPopulator {
 public:
  using SettingsCollection = Utils::UniversalSettings::DescriptorCollection;

  //! Populate Settings in the by-reference argument with default settings common to all LCAO methods.
  static void populateLcaoSettings(SettingsCollection& settings);
  //! Populate Settings in the by-reference argument with default settings common to all SCF methods.
  static void populateScfSettings(SettingsCollection& settings);
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
  static Utils::scf_mixer_t stringToScfMixer(const std::string& scfMixerName);
  /**
   * @brief This function converts a mixer type into a mixer name.
   * @param scfMixerType Type of the scf mixer.
   * @return a std::string containing the name of the scf mixer.
   */
  static std::string scfMixerToString(Utils::scf_mixer_t mixerType);

 private:
  // Lcao Settings
  static void addMolecularCharge(SettingsCollection& settings);
  static void addSpinMultiplicity(SettingsCollection& settings);
  static void addUnrestrictedCalculation(SettingsCollection& settings);
  static void addTemperatureOption(SettingsCollection& settings);
  static void addDavidsonOption(SettingsCollection& settings);

  // Scf Settings
  static void addSelfConsistanceCriterion(SettingsCollection& settings);
  static void addMaxIterations(SettingsCollection& settings);
  static void addScfMixer(SettingsCollection& settings);
};

} // namespace UniversalSettings
} // namespace Utils
} // namespace Scine

#endif // UTILS_SETTINGPOPULATOR_H
