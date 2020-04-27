/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_GAUSSIANCALCULATORSETTINGS_H
#define UTILS_GAUSSIANCALCULATORSETTINGS_H

#include <Utils/ExternalQC/Orca/OrcaCalculatorSettings.h>
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {
namespace SettingsNames {
static constexpr const char* gaussianFilenameBase = "gaussian_filename_base";
static constexpr const char* gaussianNumProcs = "gaussian_nprocs";
} // namespace SettingsNames

/**
 * @class GaussianCalculatorSettings GaussianCalculatorSettings.h
 * @brief Settings for Gaussian calculations.
 */
class GaussianCalculatorSettings : public Scine::Utils::Settings {
 public:
  // These functions populate certain settings
  void addMolecularCharge(UniversalSettings::DescriptorCollection& settings);
  void addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings);
  void addMethod(UniversalSettings::DescriptorCollection& settings);
  void addGaussianFilenameBase(UniversalSettings::DescriptorCollection& settings);
  void addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings);
  void addGaussianNumProcs(UniversalSettings::DescriptorCollection& settings);
  void addBasisSet(UniversalSettings::DescriptorCollection& settings);
  void addMemoryForGaussian(UniversalSettings::DescriptorCollection& settings);

  /**
   * @brief Constructor that populates the GaussianCalculatorSettings.
   */
  GaussianCalculatorSettings() : Settings("GaussianCalculatorSettings") {
    addMolecularCharge(_fields);
    addSpinMultiplicity(_fields);
    addMethod(_fields);
    addBasisSet(_fields);
    addGaussianFilenameBase(_fields);
    addBaseWorkingDirectory(_fields);
    addGaussianNumProcs(_fields);
    addMemoryForGaussian(_fields);
    resetToDefaults();
  };
};

inline void GaussianCalculatorSettings::addMolecularCharge(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor molecularCharge("Sets the molecular charge to use in the calculation.");
  molecularCharge.setMinimum(-10);
  molecularCharge.setMaximum(10);
  molecularCharge.setDefaultValue(0);
  settings.push_back(Scine::Utils::SettingsNames::molecularCharge, std::move(molecularCharge));
}

inline void GaussianCalculatorSettings::addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor spinMultiplicity("Sets the desired spin multiplicity to use in the "
                                                           "calculation.");
  spinMultiplicity.setMinimum(1);
  spinMultiplicity.setMaximum(10);
  spinMultiplicity.setDefaultValue(1);
  settings.push_back(Scine::Utils::SettingsNames::spinMultiplicity, std::move(spinMultiplicity));
}

inline void GaussianCalculatorSettings::addMethod(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor method("The method used in the Gaussian calculation.");
  method.setDefaultValue("PBEPBE");
  settings.push_back(Utils::SettingsNames::method, std::move(method));
}

inline void GaussianCalculatorSettings::addBasisSet(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor basisSet("The basis set used in the Gaussian calculation.");
  basisSet.setDefaultValue("def2SVP");
  settings.push_back(Utils::SettingsNames::basisSet, std::move(basisSet));
}

inline void GaussianCalculatorSettings::addGaussianFilenameBase(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor gaussianFilenameBase(
      "Base of the file name of the Gaussian calculations.");
  gaussianFilenameBase.setDefaultValue("gaussian_calc");
  settings.push_back(SettingsNames::gaussianFilenameBase, std::move(gaussianFilenameBase));
}

inline void GaussianCalculatorSettings::addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor baseWorkingDirectory("Base directory for the Gaussian calculations.");
  baseWorkingDirectory.setDefaultValue(FilesystemHelpers::currentDirectory());
  settings.push_back(SettingsNames::baseWorkingDirectory, std::move(baseWorkingDirectory));
}

inline void GaussianCalculatorSettings::addGaussianNumProcs(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor gaussianNumProcs("Number of processes for the Gaussian calculation.");
  gaussianNumProcs.setDefaultValue(1);
  gaussianNumProcs.setMinimum(1);
  settings.push_back(SettingsNames::gaussianNumProcs, std::move(gaussianNumProcs));
}

inline void GaussianCalculatorSettings::addMemoryForGaussian(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor gaussianMemory("Memory that can be used by the Gaussian calculation.");
  gaussianMemory.setDefaultValue(1024);
  settings.push_back(SettingsNames::externalQCMemory, std::move(gaussianMemory));
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_GAUSSIANCALCULATORSETTINGS_H
