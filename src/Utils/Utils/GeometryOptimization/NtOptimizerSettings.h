/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_NTOPTIMIZERSETTINGS_H_
#define UTILS_NTOPTIMIZERSETTINGS_H_

#include "Utils/GeometryOptimization/NtOptimizer.h"
#include "Utils/Settings.h"

namespace Scine {
namespace Utils {

/**
 * @brief Settings for an NtOptimizer.
 *
 * Uses template arguments in order to automatically include the
 * settings of underlying objects into the given settings.
 */
class NtOptimizerSettings : public Settings {
 public:
  /**
   * @brief Construct a new NtOptimizerSettings object.
   *
   * Sets the default values of the settings to the current values set in the objects
   * given to the constructor.
   *
   * @param nt The NT optimizer.
   */
  NtOptimizerSettings(const NtOptimizer& nt) : Settings("NtOptimizerSettings") {
    UniversalSettings::DoubleDescriptor factor("The steepest descent scaling factor.");
    factor.setDefaultValue(nt.optimizer.factor);
    factor.setMinimum(1e-12);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::sdFactor, factor);

    UniversalSettings::DoubleDescriptor repulsiveStop(
        "The stop parameter given in multiples/fractions of covalent radii sums in the repulsive case.");
    repulsiveStop.setDefaultValue(nt.check.repulsiveStop);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::repulsiveStop, repulsiveStop);

    UniversalSettings::DoubleDescriptor attractiveStop(
        "The stop parameter given in multiples/fractions of covalent radii sums in the attractive case.");
    attractiveStop.setDefaultValue(nt.check.attractiveStop);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::attractiveStop, attractiveStop);

    UniversalSettings::IntDescriptor max_iter("The maximum number of iterations.");
    max_iter.setDefaultValue(nt.check.maxIter);
    max_iter.setMinimum(1);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::maxIter, max_iter);

    UniversalSettings::IntListDescriptor nt_rhs_list(
        "The list indices of atoms to be artificially forced onto or away from those in the LHS list.");
    nt_rhs_list.setDefaultValue(nt.rhsList);
    nt_rhs_list.setItemMinimum(0);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::rHSList, nt_rhs_list);

    UniversalSettings::IntListDescriptor nt_lhs_list(
        "The list indices of atoms to be artificially forced onto or away from those in the RHS list.");
    nt_lhs_list.setDefaultValue(nt.lhsList);
    nt_lhs_list.setItemMinimum(0);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::lHSList, nt_lhs_list);

    UniversalSettings::BoolDescriptor nt_attractive("Switch for the artificial force to be attractive or repulsive.");
    nt_attractive.setDefaultValue(nt.attractive);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::attractive, nt_attractive);

    UniversalSettings::DoubleDescriptor nt_total_force_norm(
        "The norm of the summed additional forces acting on all listed atoms.");
    nt_total_force_norm.setDefaultValue(nt.totalForceNorm);
    nt_total_force_norm.setMinimum(1e-12);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::totalForceNorm, nt_total_force_norm);

    UniversalSettings::BoolDescriptor nt_use_micro_cycles(
        "Use a BFGS/GDIIS in between NT steps to run some constrained geometry optimizations.");
    nt_use_micro_cycles.setDefaultValue(nt.useMicroCycles);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::useMicroCycles, nt_use_micro_cycles);

    UniversalSettings::BoolDescriptor nt_fixed_number_of_micro_cycles(
        "Uses `numberOfMicroCycles` or grow number of micro cycles as the number of NT steps grow.");
    nt_fixed_number_of_micro_cycles.setDefaultValue(nt.fixedNumberOfMicroCycles);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::fixedNumberOfMicroCycles, nt_fixed_number_of_micro_cycles);

    UniversalSettings::IntDescriptor nt_number_of_micro_cycles("The fixed number of micro cycles.");
    nt_number_of_micro_cycles.setDefaultValue(nt.numberOfMicroCycles);
    nt_number_of_micro_cycles.setMinimum(0);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::numberOfMicroCycles, nt_number_of_micro_cycles);

    UniversalSettings::IntDescriptor nt_filter_passes(
        "Number of passes through a Savitzky-Golay filter before analyzing the reaction curve.");
    nt_filter_passes.setDefaultValue(nt.filterPasses);
    nt_filter_passes.setMinimum(0);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::filterPasses, nt_filter_passes);

    UniversalSettings::OptionListDescriptor nt_extraction("Sets the TS guess extraction criterion.");
    for (const auto& criterion : nt.possibleExtractionOptions) {
      nt_extraction.addOption(criterion);
    }
    nt_extraction.setDefaultOption(nt.possibleExtractionOptions.front());
    this->_fields.push_back(SettingsNames::Optimizations::Nt::extractionCriterion, nt_extraction);

    UniversalSettings::OptionListDescriptor nt_coordinate_system("Set the coordinate system.");
    nt_coordinate_system.addOption("internal");
    nt_coordinate_system.addOption("cartesianWithoutRotTrans");
    nt_coordinate_system.addOption("cartesian");
    nt_coordinate_system.setDefaultOption(CoordinateSystemInterpreter::getStringFromCoordinateSystem(nt.coordinateSystem));
    this->_fields.push_back(SettingsNames::Optimizations::Nt::coordinateSystem, nt_coordinate_system);

    UniversalSettings::IntListDescriptor nt_fixed_atoms("List of atoms with Cartesian constraints applied to them.");
    nt_fixed_atoms.setItemMinimum(0);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::fixedAtoms, std::move(nt_fixed_atoms));

    UniversalSettings::OptionListDescriptor nt_movable_side("Sets the sides that shall be moved.");
    nt_movable_side.addOption("lhs");
    nt_movable_side.addOption("rhs");
    nt_movable_side.addOption("both");
    nt_movable_side.setDefaultOption(nt.movableSide);
    this->_fields.push_back(SettingsNames::Optimizations::Nt::movableSide, std::move(nt_movable_side));

    this->resetToDefaults();
  }
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_NTOPTIMIZERSETTINGS_H_
