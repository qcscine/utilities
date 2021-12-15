/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Pybind.h>
#include <Utils/UniversalSettings/SettingsNames.h>

using namespace Scine::Utils;

void init_settings_names(pybind11::module& m) {
  pybind11::module settingsNames = m.def_submodule("settings_names");
  settingsNames.doc() = R"(The ``settings_names`` submodule defines common setting names to be universal across "
                         "different programs.)";

  settingsNames.attr("method") = SettingsNames::method;

  // Model
  settingsNames.attr("method_family") = SettingsNames::methodFamily;
  settingsNames.attr("spin_mode") = SettingsNames::spinMode;
  settingsNames.attr("program") = SettingsNames::program;
  settingsNames.attr("version") = SettingsNames::version;
  settingsNames.attr("basis_set") = SettingsNames::basisSet;
  settingsNames.attr("temperature") = SettingsNames::temperature;
  settingsNames.attr("electronic_temperature") = SettingsNames::electronicTemperature;
  settingsNames.attr("solvation") = SettingsNames::solvation;
  settingsNames.attr("solvent") = SettingsNames::solvent;
  settingsNames.attr("embedding") = SettingsNames::embedding;
  settingsNames.attr("periodic_boundaries") = SettingsNames::periodicBoundaries;
  settingsNames.attr("external_field") = SettingsNames::externalField;

  // General settings
  settingsNames.attr("logger_verbosity") = SettingsNames::loggerVerbosity;
  settingsNames.attr("symmetry_number") = SettingsNames::symmetryNumber;
  settingsNames.attr("method_parameters") = SettingsNames::methodParameters;
  settingsNames.attr("external_program_memory") = SettingsNames::externalProgramMemory;
  settingsNames.attr("external_program_nprocs") = SettingsNames::externalProgramNProcs;

  // SCF settings
  settingsNames.attr("molecular_charge") = SettingsNames::molecularCharge;
  settingsNames.attr("spin_multiplicity") = SettingsNames::spinMultiplicity;
  settingsNames.attr("max_scf_iterations") = SettingsNames::maxScfIterations;
  settingsNames.attr("self_consistence_criterion") = SettingsNames::selfConsistenceCriterion;
  settingsNames.attr("scf_damping") = SettingsNames::scfDamping;
  settingsNames.attr("mixer") = SettingsNames::mixer;
}
