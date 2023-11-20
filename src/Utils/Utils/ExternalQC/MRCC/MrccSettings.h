/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_MRCCSETTINGS_H
#define UTILSOS_MRCCSETTINGS_H

#include <Utils/Settings.h>

namespace Scine::Utils::ExternalQC {

/**
 * @class MrccSettings MrccSettings.h
 * @brief Settings for MRCC calculations.
 */
class MrccSettings : public Scine::Utils::Settings {
 public:
  /**
   * Constructor.
   */
  MrccSettings();

 private:
  // Populate settings.
  void addMolecularCharge(UniversalSettings::DescriptorCollection& settings);
  void addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings);
  void addSelfConsistenceCriterion(UniversalSettings::DescriptorCollection& settings);
  void addMaxScfIterations(UniversalSettings::DescriptorCollection& settings);
  void addMethod(UniversalSettings::DescriptorCollection& settings);
  void addBasisSet(UniversalSettings::DescriptorCollection& settings);
  void addSpinMode(UniversalSettings::DescriptorCollection& settings);
  void addNumProcs(UniversalSettings::DescriptorCollection& settings);
  void addMemory(UniversalSettings::DescriptorCollection& settings);
  void addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings);
  void addScfDamping(UniversalSettings::DescriptorCollection& settings);
  void addScfDampingValue(UniversalSettings::DescriptorCollection& settings);
  void addScfOrbitalShift(UniversalSettings::DescriptorCollection& settings);
  void addSolvent(UniversalSettings::DescriptorCollection& settings);
  void addSolvation(UniversalSettings::DescriptorCollection& settings);
  void addTemperature(UniversalSettings::DescriptorCollection& settings);
  void addElectronicTemperature(UniversalSettings::DescriptorCollection& settings);
  void addPressure(UniversalSettings::DescriptorCollection& settings);
};

} // namespace Scine::Utils::ExternalQC

#endif // UTILSOS_MRCCSETTINGS_H
