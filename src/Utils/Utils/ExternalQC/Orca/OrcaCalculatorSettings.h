/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ORCACALCULATORSETTINGS_H
#define UTILS_ORCACALCULATORSETTINGS_H

#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {
namespace SettingsNames {
static constexpr const char* orcaMethod = "orca_method";
static constexpr const char* orcaNumProcs = "orca_nprocs";
static constexpr const char* orcaBinaryPath = "orca_binary_path";
static constexpr const char* orcaFilenameBase = "orca_filename_base";
static constexpr const char* baseWorkingDirectory = "base_working_directory";
} // namespace SettingsNames

/**
 * @class OrcaCalculatorSettings OrcaCalculatorSettings.h
 * @brief Settings for ORCA calculations.
 */
class OrcaCalculatorSettings : public Scine::Utils::Settings {
 public:
  // These functions populate certain settings
  void addMolecularCharge(UniversalSettings::DescriptorCollection& settings);
  void addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings);
  void addSelfConsistenceCriterion(UniversalSettings::DescriptorCollection& settings);
  void addMaxIterations(UniversalSettings::DescriptorCollection& settings);
  void addOrcaMethod(UniversalSettings::DescriptorCollection& settings);
  void addOrcaNumProcs(UniversalSettings::DescriptorCollection& settings);
  void addOrcaBinaryPath(UniversalSettings::DescriptorCollection& settings);
  void addOrcaFilenameBase(UniversalSettings::DescriptorCollection& settings);
  void addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings);

  /**
   * @brief Constructor that populates the OrcaCalculatorSettings.
   */
  OrcaCalculatorSettings() : Settings("OrcaCalculatorSettings") {
    addMolecularCharge(_fields);
    addSpinMultiplicity(_fields);
    addSelfConsistenceCriterion(_fields);
    addMaxIterations(_fields);
    addOrcaMethod(_fields);
    addOrcaNumProcs(_fields);
    addOrcaBinaryPath(_fields);
    addOrcaFilenameBase(_fields);
    addBaseWorkingDirectory(_fields);
    resetToDefaults();
  };
};

inline void OrcaCalculatorSettings::addMolecularCharge(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor molecularCharge("Sets the molecular charge to use in the calculation.");
  molecularCharge.setMinimum(-10);
  molecularCharge.setMaximum(10);
  molecularCharge.setDefaultValue(0);
  settings.push_back(SettingsNames::molecularCharge, std::move(molecularCharge));
}

inline void OrcaCalculatorSettings::addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor spinMultiplicity("Sets the desired spin multiplicity to use in the "
                                                           "calculation.");
  spinMultiplicity.setMinimum(1);
  spinMultiplicity.setMaximum(10);
  spinMultiplicity.setDefaultValue(1);
  settings.push_back(SettingsNames::spinMultiplicity, std::move(spinMultiplicity));
}

inline void OrcaCalculatorSettings::addSelfConsistenceCriterion(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor selfConsistanceCriterion("Sets the desired convergence criterion.");
  selfConsistanceCriterion.setMinimum(0);
  selfConsistanceCriterion.setDefaultValue(1e-6);
  settings.push_back(SettingsNames::selfConsistanceCriterion, std::move(selfConsistanceCriterion));
}

inline void OrcaCalculatorSettings::addOrcaMethod(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor orcaMethod("The method used in the ORCA calculation.");
  orcaMethod.setDefaultValue("PBE def2-SVP");
  settings.push_back(SettingsNames::orcaMethod, std::move(orcaMethod));
}

inline void OrcaCalculatorSettings::addOrcaNumProcs(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor orcaNumProcs("Number of processes for the orca calculation.");
  orcaNumProcs.setDefaultValue(1);
  orcaNumProcs.setMinimum(1);
  settings.push_back(SettingsNames::orcaNumProcs, std::move(orcaNumProcs));
}

inline void OrcaCalculatorSettings::addOrcaBinaryPath(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor orcaBinaryPath("Path to the orca binary.");
  orcaBinaryPath.setDefaultValue("orca");
  settings.push_back(SettingsNames::orcaBinaryPath, std::move(orcaBinaryPath));
}

inline void OrcaCalculatorSettings::addOrcaFilenameBase(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor orcaFilenameBase("Base of the file name of the ORCA calculations.");
  orcaFilenameBase.setDefaultValue("orca_calc");
  settings.push_back(SettingsNames::orcaFilenameBase, std::move(orcaFilenameBase));
}

inline void OrcaCalculatorSettings::addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor baseWorkingDirectory("Base directory for the ORCA calculations.");
  baseWorkingDirectory.setDefaultValue(FilesystemHelpers::currentDirectory());
  settings.push_back(SettingsNames::baseWorkingDirectory, std::move(baseWorkingDirectory));
}

inline void OrcaCalculatorSettings::addMaxIterations(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor maxSCFIterations("Maximum number of SCF iterations.");
  maxSCFIterations.setMinimum(1);
  maxSCFIterations.setDefaultValue(100);
  settings.push_back(SettingsNames::maxIterations, std::move(maxSCFIterations));
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_ORCACALCULATORSETTINGS_H
