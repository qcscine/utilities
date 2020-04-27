/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MOLECULARDYNAMICSSETTINGS_H
#define UTILS_MOLECULARDYNAMICSSETTINGS_H

#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {

namespace SettingsNames {
static constexpr const char* timeStepInFemtoseconds = "MD_time_step";
static constexpr const char* integrationAlgorithm = "MD_integration_scheme";
static constexpr const char* temperatureBath = "temperature_bath";
static constexpr const char* targetTemperature = "target_temperature";
static constexpr const char* relaxationTimeFactor = "relaxation_time_factor";
static constexpr const char* numberOfMDSteps = "number_MD_steps";
static constexpr const char* recordFrequency = "record_frequency";
static constexpr const char* linearMomentumRemovalFrequency = "linear_momentum_removal_frequency";
static constexpr const char* angularMomentumRemovalFrequency = "angular_momentum_removal_frequency";
static constexpr const char* saveVelocities = "save_velocities";
} // namespace SettingsNames

namespace OptionNames {
static constexpr const char* leapFrogOption = "leap_frog";
static constexpr const char* eulerOption = "euler";
static constexpr const char* velocityVerletOption = "velocity_verlet";
} // namespace OptionNames

/**
 * @class MolecularDynamicsSettings MolecularDynamicsSettings.h
 * @brief Settings for Molecular Dynamics simulations.
 */
class MolecularDynamicsSettings : public Scine::Utils::Settings {
 public:
  // These functions populate certain settings
  void addTimeStep(Utils::UniversalSettings::DescriptorCollection& settings);
  void addIntegrationAlgorithm(Utils::UniversalSettings::DescriptorCollection& settings);
  void addTemperatureSettings(Utils::UniversalSettings::DescriptorCollection& settings);
  void addNumberMDSteps(Utils::UniversalSettings::DescriptorCollection& settings);
  void addRecordFrequency(Utils::UniversalSettings::DescriptorCollection& settings);
  void addMomentumRemovalFrequencies(Utils::UniversalSettings::DescriptorCollection& settings);
  void addSaveVelocitiesOption(Utils::UniversalSettings::DescriptorCollection& settings);
  /**
   * @brief Constructor that populates the MolecularDynamicsSettings.
   */
  MolecularDynamicsSettings() : Settings("MolecularDynamicsSettings") {
    addTimeStep(_fields);
    addIntegrationAlgorithm(_fields);
    addTemperatureSettings(_fields);
    addNumberMDSteps(_fields);
    addRecordFrequency(_fields);
    addMomentumRemovalFrequencies(_fields);
    addSaveVelocitiesOption(_fields);
    resetToDefaults();
  };
};

inline void MolecularDynamicsSettings::addTimeStep(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor timeStep("The MD integration time step in femtoseconds.");
  timeStep.setDefaultValue(1.0);
  settings.push_back(SettingsNames::timeStepInFemtoseconds, std::move(timeStep));
}

inline void MolecularDynamicsSettings::addIntegrationAlgorithm(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor integrationAlgorithm(
      "The integration algorithm used in the MD simulation.");
  integrationAlgorithm.addOption(OptionNames::leapFrogOption);
  integrationAlgorithm.addOption(OptionNames::eulerOption);
  integrationAlgorithm.addOption(OptionNames::velocityVerletOption);
  integrationAlgorithm.setDefaultOption(OptionNames::leapFrogOption);
  settings.push_back(SettingsNames::integrationAlgorithm, std::move(integrationAlgorithm));
}

inline void MolecularDynamicsSettings::addTemperatureSettings(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor temperatureBath(
      "Sets the coupling to a temperature bath in an MD simulation.");
  temperatureBath.setDefaultValue(true);
  settings.push_back(SettingsNames::temperatureBath, std::move(temperatureBath));

  Utils::UniversalSettings::DoubleDescriptor targetTemperature(
      "Target temperature for an MD simulation. This is only an active setting if the temperature bath coupling has "
      "been switched on.");
  targetTemperature.setMinimum(0);
  targetTemperature.setDefaultValue(300);
  settings.push_back(SettingsNames::targetTemperature, std::move(targetTemperature));

  Utils::UniversalSettings::DoubleDescriptor relaxationTimeFactor(
      "The temperature relaxation time in units of the chosen time step");
  relaxationTimeFactor.setMinimum(0.0);
  relaxationTimeFactor.setDefaultValue(10.0);
  settings.push_back(SettingsNames::relaxationTimeFactor, std::move(relaxationTimeFactor));
}

inline void MolecularDynamicsSettings::addNumberMDSteps(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor numberSteps("Number of steps in the MD simulation.");
  numberSteps.setMinimum(0);
  numberSteps.setDefaultValue(10);
  settings.push_back(SettingsNames::numberOfMDSteps, std::move(numberSteps));
}

inline void MolecularDynamicsSettings::addRecordFrequency(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor freq(
      "The frequency with which structures are written to the molecular trajectory during an MD simulation.");
  freq.setMinimum(1);
  freq.setDefaultValue(1);
  settings.push_back(SettingsNames::recordFrequency, std::move(freq));
}

inline void MolecularDynamicsSettings::addMomentumRemovalFrequencies(UniversalSettings::DescriptorCollection& settings) {
  UniversalSettings::IntDescriptor linearFreq(
      "The frequency with which the linear momentum of the center of mass is removed. If zero, no action is taken.");
  linearFreq.setMinimum(0);
  linearFreq.setDefaultValue(1);
  settings.push_back(SettingsNames::linearMomentumRemovalFrequency, std::move(linearFreq));
  UniversalSettings::IntDescriptor angularFreq(
      "The frequency with which the angular momentum of the center of mass is removed. If zero, no action is taken.");
  angularFreq.setMinimum(0);
  angularFreq.setDefaultValue(1);
  settings.push_back(SettingsNames::angularMomentumRemovalFrequency, std::move(angularFreq));
}

inline void MolecularDynamicsSettings::addSaveVelocitiesOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor saveVelocities(
      "Decides whether the velocities are saved during the MD simulation.");
  saveVelocities.setDefaultValue(false);
  settings.push_back(SettingsNames::saveVelocities, std::move(saveVelocities));
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_MOLECULARDYNAMICSSETTINGS_H
