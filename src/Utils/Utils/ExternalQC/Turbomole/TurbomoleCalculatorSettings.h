/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_TURBOMOLECALCULATORSETTINGS_H
#define UTILS_TURBOMOLECALCULATORSETTINGS_H

#include <Utils/ExternalQC/SettingsNames.h>
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @class TurbomoleCalculatorSettings TurbomoleCalculatorSettings.h
 * @brief Settings for Turbomole calculations.
 */
class TurbomoleCalculatorSettings : public Scine::Utils::Settings {
 public:
  // These functions populate certain settings
  void addMolecularCharge(UniversalSettings::DescriptorCollection& settings);
  void addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings);
  void addSelfConsistenceCriterion(UniversalSettings::DescriptorCollection& settings);
  void addMaxScfIterations(UniversalSettings::DescriptorCollection& settings);
  void addMethod(UniversalSettings::DescriptorCollection& settings);
  void addBasisSet(UniversalSettings::DescriptorCollection& settings);
  void addSpinMode(UniversalSettings::DescriptorCollection& settings);
  void addNumProcs(UniversalSettings::DescriptorCollection& settings);
  void addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings);
  void addTemperature(UniversalSettings::DescriptorCollection& settings);
  void addPressure(UniversalSettings::DescriptorCollection& settings);
  void addHessianCalculationType(UniversalSettings::DescriptorCollection& settings);
  void addElectronicTemperature(UniversalSettings::DescriptorCollection& settings);
  void addScfDamping(UniversalSettings::DescriptorCollection& settings);
  void addScfDampingValue(UniversalSettings::DescriptorCollection& settings);
  void addScfOrbitalShift(UniversalSettings::DescriptorCollection& settings);
  void addSolvent(UniversalSettings::DescriptorCollection& settings);
  void addSolvation(UniversalSettings::DescriptorCollection& settings);
  void addSteerOrbitals(UniversalSettings::DescriptorCollection& settings);
  void addPointChargesFile(UniversalSettings::DescriptorCollection& settings);
  void addEnableRi(UniversalSettings::DescriptorCollection& settings);
  void addNumExcitedStates(UniversalSettings::DescriptorCollection& settings);
  void addScfCriterionEnforce(UniversalSettings::DescriptorCollection& settings);
  void addGridEnforce(UniversalSettings::DescriptorCollection& settings);
  void addDftGrid(UniversalSettings::DescriptorCollection& settings);
  void addCavityPointsPerAtom(UniversalSettings::DescriptorCollection& settings);
  void addCavitySegmentsPerAtom(UniversalSettings::DescriptorCollection& settings);
  void addEnforceNumforce(UniversalSettings::DescriptorCollection& settings);

  /**
   * @brief Constructor that populates the TurbomoleCalculatorSettings.
   */
  TurbomoleCalculatorSettings() : Settings("TurbomoleCalculatorSettings") {
    addMolecularCharge(_fields);
    addSpinMultiplicity(_fields);
    addSelfConsistenceCriterion(_fields);
    addMaxScfIterations(_fields);
    addMethod(_fields);
    addBasisSet(_fields);
    addSpinMode(_fields);
    addNumProcs(_fields);
    addBaseWorkingDirectory(_fields);
    addTemperature(_fields);
    addPressure(_fields);
    addScfDamping(_fields);
    addScfDampingValue(_fields);
    addScfOrbitalShift(_fields);
    addHessianCalculationType(_fields);
    addElectronicTemperature(_fields);
    addSolvent(_fields);
    addSolvation(_fields);
    addSteerOrbitals(_fields);
    addPointChargesFile(_fields);
    addEnableRi(_fields);
    addNumExcitedStates(_fields);
    addScfCriterionEnforce(_fields);
    addGridEnforce(_fields);
    addDftGrid(_fields);
    addCavityPointsPerAtom(_fields);
    addCavitySegmentsPerAtom(_fields);
    addEnforceNumforce(_fields);
    resetToDefaults();
  };
};

inline void TurbomoleCalculatorSettings::addMolecularCharge(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor molecularCharge("Sets the molecular charge to use in the calculation.");
  molecularCharge.setMinimum(-10);
  molecularCharge.setMaximum(10);
  molecularCharge.setDefaultValue(0);
  settings.push_back(Scine::Utils::SettingsNames::molecularCharge, std::move(molecularCharge));
}

inline void TurbomoleCalculatorSettings::addSpinMultiplicity(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor spinMultiplicity("Sets the desired spin multiplicity to use in the "
                                                           "calculation.");
  spinMultiplicity.setMinimum(1);
  spinMultiplicity.setMaximum(10);
  spinMultiplicity.setDefaultValue(1);
  settings.push_back(Scine::Utils::SettingsNames::spinMultiplicity, std::move(spinMultiplicity));
}

inline void TurbomoleCalculatorSettings::addSelfConsistenceCriterion(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor selfConsistenceCriterion("Sets the desired convergence criterion.");
  selfConsistenceCriterion.setMinimum(0);
  selfConsistenceCriterion.setDefaultValue(1e-7);
  settings.push_back(Scine::Utils::SettingsNames::selfConsistenceCriterion, std::move(selfConsistenceCriterion));
}

inline void TurbomoleCalculatorSettings::addMethod(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor method("The method used in the Turbomole calculation.");
  method.setDefaultValue("pbe");
  settings.push_back(Utils::SettingsNames::method, std::move(method));
}

inline void TurbomoleCalculatorSettings::addBasisSet(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor basisSet("The basis set used in the Turbomole calculation.");
  basisSet.setDefaultValue("def2-SVP");
  settings.push_back(Utils::SettingsNames::basisSet, std::move(basisSet));
}

inline void TurbomoleCalculatorSettings::addSpinMode(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor spinMode("The spin mode such as 'restricted' or 'unrestricted'.");
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Any));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Restricted));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::RestrictedOpenShell));
  spinMode.addOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Unrestricted));
  spinMode.setDefaultOption(SpinModeInterpreter::getStringFromSpinMode(SpinMode::Any));
  settings.push_back(Utils::SettingsNames::spinMode, std::move(spinMode));
}

inline void TurbomoleCalculatorSettings::addNumProcs(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor turbomoleNumProcs("Number of processes for the Turbomole calculation.");
  turbomoleNumProcs.setDefaultValue(1);
  turbomoleNumProcs.setMinimum(1);
  settings.push_back(Utils::SettingsNames::externalProgramNProcs, std::move(turbomoleNumProcs));
}

inline void TurbomoleCalculatorSettings::addBaseWorkingDirectory(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor baseWorkingDirectory("Base directory for the Turbomole calculations.");
  baseWorkingDirectory.setDefaultValue(FilesystemHelpers::currentDirectory());
  settings.push_back(SettingsNames::baseWorkingDirectory, std::move(baseWorkingDirectory));
}

inline void TurbomoleCalculatorSettings::addMaxScfIterations(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor maxScfIterations("Maximum number of SCF iterations.");
  maxScfIterations.setMinimum(1);
  maxScfIterations.setDefaultValue(100);
  settings.push_back(Scine::Utils::SettingsNames::maxScfIterations, std::move(maxScfIterations));
}

inline void TurbomoleCalculatorSettings::addTemperature(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor temperature("Sets the temperature for the thermochemical calculation.");
  temperature.setDefaultValue(298.15);
  settings.push_back(Utils::SettingsNames::temperature, std::move(temperature));
}

inline void TurbomoleCalculatorSettings::addPressure(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor pressure("Sets the pressure for the thermochemical calculation in Pa.");
  pressure.setDefaultValue(101325.0);
  settings.push_back(Utils::SettingsNames::pressure, std::move(pressure));
}

inline void TurbomoleCalculatorSettings::addScfDamping(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor scfDamping("Enable stronger SCF damping (true/false).");
  scfDamping.setDefaultValue(false);
  settings.push_back(Utils::SettingsNames::scfDamping, std::move(scfDamping));
}

inline void TurbomoleCalculatorSettings::addScfDampingValue(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor scfDampingValue("Specify exact SCF damping value to be used.");
  scfDampingValue.setDefaultValue(0.5);
  settings.push_back(SettingsNames::scfDampingValue, std::move(scfDampingValue));
}

inline void TurbomoleCalculatorSettings::addScfOrbitalShift(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor scfOrbitalShift(
      "Shift virtual orbitals to higher energies to aid convergence.");
  // internal default
  scfOrbitalShift.setDefaultValue(0.1);
  settings.push_back(SettingsNames::scfOrbitalShift, std::move(scfOrbitalShift));
}

inline void TurbomoleCalculatorSettings::addHessianCalculationType(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor hessianType("The method for calculating the Hessian.");
  hessianType.addOption("analytical");
  hessianType.addOption("numerical");
  hessianType.setDefaultOption("analytical");
  settings.push_back(ExternalQC::SettingsNames::hessianCalculationType, std::move(hessianType));
}

inline void TurbomoleCalculatorSettings::addElectronicTemperature(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor electronicTemperature(
      "Sets the electronic temperature for SCF calculations.");
  electronicTemperature.setMinimum(0.0);
  electronicTemperature.setDefaultValue(0.0);
  settings.push_back(Utils::SettingsNames::electronicTemperature, std::move(electronicTemperature));
}

inline void TurbomoleCalculatorSettings::addSolvent(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor solventOption("Sets the implicit solvent.");
  solventOption.setDefaultValue("");
  settings.push_back(Utils::SettingsNames::solvent, std::move(solventOption));
}

inline void TurbomoleCalculatorSettings::addSolvation(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor solvationOption(
      "Sets the implicit solvation model in the TURBOMOLE calculation.");
  solvationOption.setDefaultValue("");
  settings.push_back(Utils::SettingsNames::solvation, std::move(solvationOption));
}

inline void TurbomoleCalculatorSettings::addSteerOrbitals(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor steerOrbitals(
      "Converts internal coordinates used by default to cartesian coordinates.");
  steerOrbitals.setDefaultValue(false);
  settings.push_back(SettingsNames::steerOrbitals, std::move(steerOrbitals));
}

inline void TurbomoleCalculatorSettings::addPointChargesFile(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor pointChargesFile(
      "Sets the file name for a Turbomole point charges file. Note that the expected line format for the point "
      "charges file is <x> <y> <z> <q>.");
  pointChargesFile.setDefaultValue("");
  settings.push_back(SettingsNames::pointChargesFile, std::move(pointChargesFile));
}

inline void TurbomoleCalculatorSettings::addEnableRi(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor enableRi("Enables the Resolution of the Identity Approximation.");
  enableRi.setDefaultValue(true);
  settings.push_back(SettingsNames::enableRi, std::move(enableRi));
}

inline void TurbomoleCalculatorSettings::addNumExcitedStates(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor numExcitedStates(
      "The total number of electronically excited states to be calculated. Note that properties such as energy "
      "and nuclear gradients are only calculated for the highest excited state.");
  numExcitedStates.setDefaultValue(0);
  numExcitedStates.setMinimum(0);
  settings.push_back(SettingsNames::numExcitedStates, std::move(numExcitedStates));
}

inline void TurbomoleCalculatorSettings::addScfCriterionEnforce(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor scfEnforce(
      "Whether the set self_consistence_criterion should not be made stricter, "
      "even if derivative quantities are calculated.");
  scfEnforce.setDefaultValue(false);
  settings.push_back(SettingsNames::enforceScfCriterion, std::move(scfEnforce));
}

inline void TurbomoleCalculatorSettings::addGridEnforce(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor gridEnforce(
      "If true, the grid accuracy is used as provided in the input"
      " even for Hessian calculations. If false, the grid accuracy is increased to"
      " m4 if the grid is considered too coarse for Hessian calculations.");
  gridEnforce.setDefaultValue(false);
  settings.push_back(SettingsNames::enforceGrid, std::move(gridEnforce));
}

inline void TurbomoleCalculatorSettings::addDftGrid(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor dftGrid(
      "Specify DFT grid to be used."
      "Possible grids range from 1-7 and m3-m5, respectively, where 1 is coarse and 7 most dense.");
  // Add all possible options
  dftGrid.addOption("m3");
  dftGrid.addOption("m4");
  dftGrid.addOption("m5");
  for (int i = 1; i < 8; i++) {
    dftGrid.addOption(std::to_string(i));
  }
  // Set Default value
  dftGrid.setDefaultOption("m3");
  settings.push_back(SettingsNames::dftGrid, std::move(dftGrid));
}

inline void TurbomoleCalculatorSettings::addCavityPointsPerAtom(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor cavityPointsPerAtom(
      "The number basis grid points per atom for the cavity construction"
      "Allowed values must fulfill: i = 10 * 3^k * 4^l + 2");
  cavityPointsPerAtom.setDefaultValue(1082);
  cavityPointsPerAtom.setMinimum(12);
  settings.push_back(SettingsNames::cavityPointsPerAtom, std::move(cavityPointsPerAtom));
}

inline void TurbomoleCalculatorSettings::addCavitySegmentsPerAtom(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor cavitySegmentsPerAtom(
      "The number of segments per atom for the cavity construction"
      "Allowed values must fulfill: i = 10 * 3^k * 4^l + 2");
  cavitySegmentsPerAtom.setDefaultValue(92);
  cavitySegmentsPerAtom.setMinimum(12);
  settings.push_back(SettingsNames::cavitySegmentsPerAtom, std::move(cavitySegmentsPerAtom));
}

inline void TurbomoleCalculatorSettings::addEnforceNumforce(UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor enforceNumforce(
      "Whether Turbomole should skip its gradient check when performing numforce.");
  enforceNumforce.setDefaultValue(false);
  settings.push_back(SettingsNames::enforceNumforce, std::move(enforceNumforce));
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_TURBOMOLECALCULATORSETTINGS_H
