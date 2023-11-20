/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_AFIROPTIMIZERSETTINGS_H_
#define UTILS_AFIROPTIMIZERSETTINGS_H_

#include "Utils/GeometryOptimization/AfirOptimizerBase.h"
#include "Utils/Settings.h"
#include "Utils/UniversalSettings/OptimizationSettingsNames.h"

namespace Scine {
namespace Utils {

/**
 * @brief Settings for an AfirOptimizer.
 *
 * Uses template arguments in order to automatically include the
 * settings of underlying objects into the given settings.
 *
 * @tparam OptimizerType The underlying Optimizer class.
 * @tparam ConvergenceCheckType The underlying ConvergenceCheck class.
 */
template<class OptimizerType, class ConvergenceCheckType>
class AfirOptimizerSettings : public Settings {
 public:
  /**
   * @brief Construct a new AfirOptimizerSettings object.
   *
   * Sets the default values of the settings to the current values set in the objects
   * given to the constructor.
   *
   * @param afirBase The AFIR optimizer.
   * @param optimizer The optimizer.
   * @param check The convergence check criteria.
   */
  AfirOptimizerSettings(const AfirOptimizerBase& afirBase, const OptimizerType& optimizer, const ConvergenceCheckType& check)
    : Settings("AfirOptimizerSettings") {
    optimizer.addSettingsDescriptors(this->_fields);
    check.addSettingsDescriptors(this->_fields);
    check.addAfirSettingsDescriptors(this->_fields);

    UniversalSettings::IntListDescriptor afir_rhs_list(
        "The list of atoms indices of atoms to be artificially forced onto or away from those in the LHS list.");
    afir_rhs_list.setDefaultValue(afirBase.rhsList);
    this->_fields.push_back(SettingsNames::Optimizations::Afir::rHSList, afir_rhs_list);

    UniversalSettings::IntListDescriptor afir_lhs_list(
        "The list of atoms indices of atoms to be artificially forced onto or away from those in the RHS list.");
    afir_lhs_list.setDefaultValue(afirBase.lhsList);
    this->_fields.push_back(SettingsNames::Optimizations::Afir::lHSList, afir_lhs_list);

    UniversalSettings::BoolDescriptor afir_weak_forces(
        "Switch for an additional weak attractive force applied to all atom pairs.");
    afir_weak_forces.setDefaultValue(afirBase.weak);
    this->_fields.push_back(SettingsNames::Optimizations::Afir::weakForces, afir_weak_forces);

    UniversalSettings::BoolDescriptor afir_attractive("Switch for the artificial force to be attractive or repulsive.");
    afir_attractive.setDefaultValue(afirBase.attractive);
    this->_fields.push_back(SettingsNames::Optimizations::Afir::attractive, afir_attractive);

    UniversalSettings::DoubleDescriptor afir_energy_allowance(
        "The maximum amount of energy to be added by the artifical force, in kJ/mol.");
    afir_energy_allowance.setDefaultValue(afirBase.energyAllowance);
    this->_fields.push_back(SettingsNames::Optimizations::Afir::energyAllowance, afir_energy_allowance);

    UniversalSettings::IntDescriptor afir_phase_in(
        "The number of steps over which the full attractive force is slowly applied");
    afir_phase_in.setDefaultValue(afirBase.phaseIn);
    this->_fields.push_back(SettingsNames::Optimizations::Afir::phaseIn, afir_phase_in);

    UniversalSettings::OptionListDescriptor afir_coordinate_system("Set the coordinate system.");
    afir_coordinate_system.addOption("internal");
    afir_coordinate_system.addOption("cartesianWithoutRotTrans");
    afir_coordinate_system.addOption("cartesian");
    afir_coordinate_system.setDefaultOption(
        CoordinateSystemInterpreter::getStringFromCoordinateSystem(afirBase.coordinateSystem));
    this->_fields.push_back(SettingsNames::Optimizations::Afir::coordinateSystem, afir_coordinate_system);

    this->resetToDefaults();
  }
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_AFIROPTIMIZERSETTINGS_H_
