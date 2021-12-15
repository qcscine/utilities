/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_LENNARDJONESCALCULATORSETTINGS_H
#define UTILS_LENNARDJONESCALCULATORSETTINGS_H

#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {
namespace SettingsNames {
static constexpr const char* lennardJonesSigma = "lj_sigma";
static constexpr const char* lennardJonesEpsilon = "lj_epsilon";
static constexpr const char* lennardJonesCutoff = "lj_cutoff";
static constexpr const char* lennardJonesUsePBCs = "lj_use_pbcs";
static constexpr const char* lennardJonesBoxsize = "lj_boxsize";
} // namespace SettingsNames

/**
 * @class LennardJonesCalculatorSettings LennardJonesCalculatorSettings.h
 * @brief Settings for Lennard-Jones calculations.
 */
class LennardJonesCalculatorSettings : public Scine::Utils::Settings {
 public:
  /**
   * @brief Populates the Lennard-Jones settings
   */
  void populateSettings(UniversalSettings::DescriptorCollection& settings);

  /**
   * @brief Constructor that populates the LennardJonesCalculatorSettings.
   */
  LennardJonesCalculatorSettings() : Settings("LennardJonesCalculatorSettings") {
    populateSettings(_fields);
    resetToDefaults();
  };
};

inline void LennardJonesCalculatorSettings::populateSettings(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor lennardJonesSigma(
      "The sigma parameter for the Lennard-Jones potential in bohr.");
  lennardJonesSigma.setMinimum(0.);
  lennardJonesSigma.setDefaultValue(6.4); // reasonable for Ar
  settings.push_back(Scine::Utils::SettingsNames::lennardJonesSigma, std::move(lennardJonesSigma));

  Utils::UniversalSettings::DoubleDescriptor lennardJonesEpsilon(
      "The depth epsilon of the Lennard-Jones potential in K.");
  lennardJonesEpsilon.setMinimum(0.);
  lennardJonesEpsilon.setDefaultValue(120); // reasonable for Ar
  settings.push_back(Scine::Utils::SettingsNames::lennardJonesEpsilon, std::move(lennardJonesEpsilon));

  Utils::UniversalSettings::DoubleDescriptor lennardJonesCutoff(
      "The cutoff radius for the Lennard-Jones potential in bohr.");
  lennardJonesCutoff.setMinimum(0.);
  lennardJonesCutoff.setDefaultValue(16.);
  settings.push_back(Scine::Utils::SettingsNames::lennardJonesCutoff, std::move(lennardJonesCutoff));

  Utils::UniversalSettings::BoolDescriptor lennardJonesUsePBCs("Whether periodic boundary conditions are to be used.");
  lennardJonesUsePBCs.setDefaultValue(false);
  settings.push_back(Scine::Utils::SettingsNames::lennardJonesUsePBCs, std::move(lennardJonesUsePBCs));

  Utils::UniversalSettings::DoubleDescriptor lennardJonesBoxsize(
      "The side length of the cubic simulation cell in bohr. Has to be more than twice as large as cut-off radius. "
      "Only active if periodic boundary conditions are turned on.");
  lennardJonesBoxsize.setMinimum(0.);
  lennardJonesBoxsize.setDefaultValue(40);
  settings.push_back(Scine::Utils::SettingsNames::lennardJonesBoxsize, std::move(lennardJonesBoxsize));
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_LENNARDJONESCALCULATORSETTINGS_H
