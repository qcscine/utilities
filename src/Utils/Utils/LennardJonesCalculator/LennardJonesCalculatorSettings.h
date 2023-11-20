/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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
  UniversalSettings::DoubleDescriptor convergence_threshold("Energy convergence limit.");
  convergence_threshold.setDefaultValue(1e-12);
  settings.push_back(Scine::Utils::SettingsNames::selfConsistenceCriterion, convergence_threshold);

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

  Utils::UniversalSettings::StringDescriptor lennardJonesPBCs(
      "The periodic boundary conditions. Empty if not applied.");
  lennardJonesPBCs.setDefaultValue("");
  settings.push_back(Scine::Utils::SettingsNames::periodicBoundaries, std::move(lennardJonesPBCs));
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_LENNARDJONESCALCULATORSETTINGS_H
