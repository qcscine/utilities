/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_NTOPTIMIZER2SETTINGS_H_
#define UTILS_NTOPTIMIZER2SETTINGS_H_

#include "Utils/GeometryOptimization/NtOptimizer2.h"
#include "Utils/Settings.h"

namespace Scine {
namespace Utils {

/**
 * @brief Settings for an NtOptimizer2.
 *
 * Uses template arguments in order to automatically include the
 * settings of underlying objects into the given settings.
 */
class NtOptimizer2Settings : public Settings {
 public:
  /**
   * @brief Construct a new NtOptimizer2Settings object.
   *
   * Sets the default values of the settings to the current values set in the objects
   * given to the constructor.
   *
   * @param nt The NTOptimizer2.
   */
  NtOptimizer2Settings(const NtOptimizer2& nt) : Settings("NtOptimizer2Settings") {
    UniversalSettings::DoubleDescriptor factor("The steepest descent scaling factor.");
    factor.setDefaultValue(nt.optimizer.factor);
    factor.setMinimum(1e-12);
    this->_fields.push_back(NtOptimizer2::ntSdFactorKey, factor);

    UniversalSettings::DoubleDescriptor attractiveStop(
        "The stop parameter given in multiples/fractions of covalent radii sums in the attractive case.");
    attractiveStop.setDefaultValue(nt.check.attractiveDistanceStop);
    this->_fields.push_back(NtOptimizer2::ntAttractiveStopKey, attractiveStop);

    UniversalSettings::IntDescriptor max_iter("The maximum number of iterations.");
    max_iter.setDefaultValue(nt.check.maxIter);
    max_iter.setMinimum(1);
    this->_fields.push_back(NtOptimizer2::ntMaxIterKey, max_iter);

    UniversalSettings::IntListDescriptor nt_associations("List of atom pairs to be pushed together to for a bond.");
    nt_associations.setDefaultValue(nt.associationList);
    nt_associations.setItemMinimum(0);
    this->_fields.push_back(NtOptimizer2::ntAssListKey, nt_associations);

    UniversalSettings::IntListDescriptor nt_dissociations(
        "List of atom pairs to be pulled apart, breaking their bond.");
    nt_dissociations.setDefaultValue(nt.dissociationList);
    nt_dissociations.setItemMinimum(0);
    this->_fields.push_back(NtOptimizer2::ntDissListKey, nt_dissociations);

    UniversalSettings::DoubleDescriptor nt_total_force_norm(
        "The norm of the summed additional forces acting on all listed atoms.");
    nt_total_force_norm.setDefaultValue(nt.totalForceNorm);
    nt_total_force_norm.setMinimum(1e-12);
    this->_fields.push_back(NtOptimizer2::ntTotalForceNormKey, nt_total_force_norm);

    UniversalSettings::BoolDescriptor nt_use_micro_cycles(
        "Use a BFGS/GDIIS in between NT steps to run some constrained geometry optimizations.");
    nt_use_micro_cycles.setDefaultValue(nt.useMicroCycles);
    this->_fields.push_back(NtOptimizer2::ntUseMicroCycles, nt_use_micro_cycles);

    UniversalSettings::BoolDescriptor nt_fixed_number_of_micro_cycles(
        "Uses `numberOfMicroCycles` or grow number of micro cycles as the number of NT steps grow.");
    nt_fixed_number_of_micro_cycles.setDefaultValue(nt.fixedNumberOfMicroCycles);
    this->_fields.push_back(NtOptimizer2::ntFixedNumberOfMicroCycles, nt_fixed_number_of_micro_cycles);

    UniversalSettings::IntDescriptor nt_number_of_micro_cycles("The fixed number of micro cycles.");
    nt_number_of_micro_cycles.setDefaultValue(nt.numberOfMicroCycles);
    nt_number_of_micro_cycles.setMinimum(0);
    this->_fields.push_back(NtOptimizer2::ntNumberOfMicroCycles, nt_number_of_micro_cycles);

    UniversalSettings::IntDescriptor nt_filter_passes(
        "Number of passes through a Savitzky-Golay filter before analyzing the reaction curve.");
    nt_filter_passes.setDefaultValue(nt.filterPasses);
    nt_filter_passes.setMinimum(0);
    this->_fields.push_back(NtOptimizer2::ntFilterPasses, nt_filter_passes);

    UniversalSettings::OptionListDescriptor nt_coordinate_system("Set the coordinate system.");
    nt_coordinate_system.addOption("internal");
    nt_coordinate_system.addOption("cartesianWithoutRotTrans");
    nt_coordinate_system.addOption("cartesian");
    nt_coordinate_system.setDefaultOption(CoordinateSystemInterpreter::getStringFromCoordinateSystem(nt.coordinateSystem));
    this->_fields.push_back(NtOptimizer2::ntCoordinateSystemKey, nt_coordinate_system);

    UniversalSettings::IntListDescriptor nt_fixed_atoms("List of atoms with Cartesian constraints applied to them.");
    nt_fixed_atoms.setItemMinimum(0);
    this->_fields.push_back(NtOptimizer2::ntFixedAtomsKey, std::move(nt_fixed_atoms));

    this->resetToDefaults();
  }
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_NTOPTIMIZER2SETTINGS_H_
