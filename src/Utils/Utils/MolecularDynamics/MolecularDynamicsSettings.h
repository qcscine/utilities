/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MOLECULARDYNAMICSSETTINGS_H
#define UTILS_MOLECULARDYNAMICSSETTINGS_H

#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {

namespace SettingsNames {
static constexpr const char* timeStepInFemtoseconds = "md_time_step";
static constexpr const char* integrationAlgorithm = "md_integration_scheme";
static constexpr const char* generationTemperature = "generation_temperature";
static constexpr const char* generationSeed = "generation_seed";
static constexpr const char* thermostatAlgorithm = "md_thermostat";
static constexpr const char* targetTemperature = "target_temperature";
static constexpr const char* temperatureCouplingTime = "temperature_coupling_time";
static constexpr const char* stochasticDynamicsSeed = "stochastic_dynamics_seed";
static constexpr const char* numberOfMDSteps = "number_md_steps";
static constexpr const char* recordFrequency = "record_frequency";
static constexpr const char* linearMomentumRemovalFrequency = "linear_momentum_removal_frequency";
static constexpr const char* angularMomentumRemovalFrequency = "angular_momentum_removal_frequency";
static constexpr const char* saveVelocities = "save_velocities";
static constexpr const char* saveTemperatures = "save_temperatures";
static constexpr const char* requireCharges = "require_charges";
static constexpr const char* requireBondOrders = "require_bond_orders";
} // namespace SettingsNames

namespace OptionNames {
// Integrators
static constexpr const char* leapFrogOption = "leap_frog";
static constexpr const char* eulerOption = "euler";
static constexpr const char* velocityVerletOption = "velocity_verlet";
static constexpr const char* stochasticDynamicsOption = "stochastic_dynamics";
// Thermostat names
static constexpr const char* noThermostatOption = "none";
static constexpr const char* berendsenThermostatOption = "berendsen";
} // namespace OptionNames

// Default values for the temperature coupling times in femtoseconds
namespace TemperatureCouplingDefaults {
static constexpr const double berendsenTimeDefault = 10;            // Sensible for equilibration, larger for production
static constexpr const double stochasticDynamicsTimeDefault = 2000; // 1/gamma
} // namespace TemperatureCouplingDefaults

/**
 * @class MolecularDynamicsSettings MolecularDynamicsSettings.h
 * @brief Settings for Molecular Dynamics simulations.
 */
class MolecularDynamicsSettings : public Scine::Utils::Settings {
 public:
  // These functions populate certain settings
  void addGenerationSettings(Utils::UniversalSettings::DescriptorCollection& settings);
  void addTimeStep(Utils::UniversalSettings::DescriptorCollection& settings);
  void addIntegrationAlgorithm(Utils::UniversalSettings::DescriptorCollection& settings);
  void addTemperatureSettings(Utils::UniversalSettings::DescriptorCollection& settings);
  void addNumberMDSteps(Utils::UniversalSettings::DescriptorCollection& settings);
  void addRecordFrequency(Utils::UniversalSettings::DescriptorCollection& settings);
  void addMomentumRemovalFrequencies(Utils::UniversalSettings::DescriptorCollection& settings);
  void addSaveVelocitiesOption(Utils::UniversalSettings::DescriptorCollection& settings);
  void addSaveTemperaturesOption(Utils::UniversalSettings::DescriptorCollection& settings);
  void addCalculatorProperties(Utils::UniversalSettings::DescriptorCollection& settings);
  /**
   * @brief Constructor that populates the MolecularDynamicsSettings.
   */
  MolecularDynamicsSettings() : Settings("MolecularDynamicsSettings") {
    addGenerationSettings(_fields);
    addTimeStep(_fields);
    addIntegrationAlgorithm(_fields);
    addTemperatureSettings(_fields);
    addNumberMDSteps(_fields);
    addRecordFrequency(_fields);
    addMomentumRemovalFrequencies(_fields);
    addSaveVelocitiesOption(_fields);
    addSaveTemperaturesOption(_fields);
    addCalculatorProperties(_fields);
    resetToDefaults();
  };
};

// Settings for the inital velocity generation
inline void MolecularDynamicsSettings::addGenerationSettings(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor generationTemperature(
      "Temperature in K for which initial velocities are drawn from a Boltzmann distribution, unless they are given "
      "explicitly. If zero, all initial velocities are set to zero.");
  generationTemperature.setMinimum(0);
  generationTemperature.setDefaultValue(300);
  settings.push_back(Utils::SettingsNames::generationTemperature, std::move(generationTemperature));

  Utils::UniversalSettings::IntDescriptor generationSeed("The seed to draw the initial velocity distribution.");
  generationSeed.setDefaultValue(42);
  settings.push_back(Utils::SettingsNames::generationSeed, std::move(generationSeed));
}

inline void MolecularDynamicsSettings::addTimeStep(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor timeStep("The MD integration time step in femtoseconds.");
  timeStep.setDefaultValue(1.0);
  settings.push_back(Utils::SettingsNames::timeStepInFemtoseconds, std::move(timeStep));
}

inline void MolecularDynamicsSettings::addIntegrationAlgorithm(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor integrationAlgorithm(
      "The integration algorithm used in the MD simulation.");
  integrationAlgorithm.addOption(OptionNames::leapFrogOption);
  integrationAlgorithm.addOption(OptionNames::eulerOption);
  integrationAlgorithm.addOption(OptionNames::velocityVerletOption);
  integrationAlgorithm.addOption(OptionNames::stochasticDynamicsOption);
  integrationAlgorithm.setDefaultOption(OptionNames::leapFrogOption);
  settings.push_back(Utils::SettingsNames::integrationAlgorithm, std::move(integrationAlgorithm));
}

inline void MolecularDynamicsSettings::addTemperatureSettings(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor thermostatAlgorithm(
      "Sets the coupling to a temperature bath in an MD simulation.");
  thermostatAlgorithm.addOption(OptionNames::berendsenThermostatOption);
  thermostatAlgorithm.addOption(OptionNames::noThermostatOption);
  thermostatAlgorithm.setDefaultOption(OptionNames::noThermostatOption);
  settings.push_back(Utils::SettingsNames::thermostatAlgorithm, std::move(thermostatAlgorithm));

  Utils::UniversalSettings::DoubleDescriptor targetTemperature(
      "Target temperature in K for an MD simulation. If zero, the generation temperature is used."
      "This is only an active setting with stochastic dynamics or a thermostat.");
  targetTemperature.setMinimum(0);
  targetTemperature.setDefaultValue(0);
  settings.push_back(Utils::SettingsNames::targetTemperature, std::move(targetTemperature));

  Utils::UniversalSettings::DoubleDescriptor temperatureCouplingTime(
      "The thermostat time parameter in fs. If set to zero the default parameter of the chosen thermostat is used.");
  temperatureCouplingTime.setMinimum(0.0);
  temperatureCouplingTime.setDefaultValue(0.0);
  settings.push_back(Utils::SettingsNames::temperatureCouplingTime, std::move(temperatureCouplingTime));

  Utils::UniversalSettings::IntDescriptor stochasticDynamicsSeed("The seed used for stochastic dynamics.");
  stochasticDynamicsSeed.setDefaultValue(42);
  settings.push_back(Utils::SettingsNames::stochasticDynamicsSeed, std::move(stochasticDynamicsSeed));
}

inline void MolecularDynamicsSettings::addNumberMDSteps(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor numberSteps("Number of steps in the MD simulation.");
  numberSteps.setMinimum(0);
  numberSteps.setDefaultValue(10);
  settings.push_back(Utils::SettingsNames::numberOfMDSteps, std::move(numberSteps));
}

inline void MolecularDynamicsSettings::addRecordFrequency(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor freq(
      "The frequency with which structures are written to the molecular trajectory during an MD simulation.");
  freq.setMinimum(1);
  freq.setDefaultValue(1);
  settings.push_back(Utils::SettingsNames::recordFrequency, std::move(freq));
}

inline void MolecularDynamicsSettings::addMomentumRemovalFrequencies(UniversalSettings::DescriptorCollection& settings) {
  UniversalSettings::IntDescriptor linearFreq(
      "The frequency with which the linear momentum of the center of mass is removed. If zero, no action is taken.");
  linearFreq.setMinimum(0);
  linearFreq.setDefaultValue(1);
  settings.push_back(Utils::SettingsNames::linearMomentumRemovalFrequency, std::move(linearFreq));
  UniversalSettings::IntDescriptor angularFreq(
      "The frequency with which the angular momentum of the center of mass is removed. If zero, no action is taken.");
  angularFreq.setMinimum(0);
  angularFreq.setDefaultValue(1);
  settings.push_back(Utils::SettingsNames::angularMomentumRemovalFrequency, std::move(angularFreq));
}

inline void MolecularDynamicsSettings::addSaveVelocitiesOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor saveVelocities(
      "Decides whether the velocities are saved during the MD simulation.");
  saveVelocities.setDefaultValue(false);
  settings.push_back(Utils::SettingsNames::saveVelocities, std::move(saveVelocities));
}

inline void MolecularDynamicsSettings::addSaveTemperaturesOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor saveTemperatures(
      "Decides whether the temperatures are saved during the MD simulation.");
  saveTemperatures.setDefaultValue(false);
  settings.push_back(Utils::SettingsNames::saveTemperatures, std::move(saveTemperatures));
}

inline void MolecularDynamicsSettings::addCalculatorProperties(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor requireCharges(
      "Whether the calculator shall calculate charges during the MD simulation.");
  requireCharges.setDefaultValue(false);
  settings.push_back(Utils::SettingsNames::requireCharges, std::move(requireCharges));

  Utils::UniversalSettings::BoolDescriptor requireBondOrders(
      "Whether the calculator shall calculate bond orders during the MD simulation.");
  requireBondOrders.setDefaultValue(false);
  settings.push_back(Utils::SettingsNames::requireBondOrders, std::move(requireBondOrders));
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_MOLECULARDYNAMICSSETTINGS_H
