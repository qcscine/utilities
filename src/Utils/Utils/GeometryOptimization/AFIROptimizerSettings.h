/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_AFIROPTIMIZERSETTINGS_H_
#define UTILS_AFIROPTIMIZERSETTINGS_H_

#include "Utils/GeometryOptimization/AFIROptimizerBase.h"
#include "Utils/Settings.h"

namespace Scine {
namespace Utils {

/**
 * @brief Settings for an AFIROptimizer.
 *
 * Uses template arguments in order to automatically include the
 * settings of underlying objects into the given settings.
 *
 * @tparam OptimizerType The underlying Optimizer class.
 * @tparam ConvergenceCheckType The underlying ConvergenceCheck class.
 */
template<class OptimizerType, class ConvergenceCheckType>
class AFIROptimizerSettings : public Settings {
 public:
  /**
   * @brief Construct a new AFIROptimizerSettings object.
   *
   * Sets the default values of the settings to the current values set in the objects
   * given to the constructor.
   *
   * @param afirBase The AFIR optimizer.
   * @param optimizer The optimizer.
   * @param check The convergence check criteria.
   */
  AFIROptimizerSettings(const AFIROptimizerBase& afirBase, const OptimizerType& optimizer, const ConvergenceCheckType& check)
    : Settings("AFIROptimizerSettings") {
    optimizer.addSettingsDescriptors(this->_fields);
    check.addSettingsDescriptors(this->_fields);

    UniversalSettings::IntListDescriptor afir_rhs_list(
        "The list of atoms indices of atoms to be artificially forced onto or away from those in the LHS list.");
    afir_rhs_list.setDefaultValue(afirBase.rhsList);
    this->_fields.push_back(AFIROptimizerBase::afirRHSListKey, afir_rhs_list);

    UniversalSettings::IntListDescriptor afir_lhs_list(
        "The list of atoms indices of atoms to be artificially forced onto or away from those in the RHS list.");
    afir_lhs_list.setDefaultValue(afirBase.lhsList);
    this->_fields.push_back(AFIROptimizerBase::afirLHSListKey, afir_lhs_list);

    UniversalSettings::BoolDescriptor afir_weak_forces(
        "Switch for an additional weak attractive force applied to all atom pairs.");
    afir_weak_forces.setDefaultValue(afirBase.weak);
    this->_fields.push_back(AFIROptimizerBase::afirWeakForcesKey, afir_weak_forces);

    UniversalSettings::BoolDescriptor afir_attractive("Switch for the artificial force to be attractive or repulsive.");
    afir_attractive.setDefaultValue(afirBase.attractive);
    this->_fields.push_back(AFIROptimizerBase::afirAttractiveKey, afir_attractive);

    UniversalSettings::DoubleDescriptor afir_energy_allowance(
        "The maximum amount of energy to be added by the artifical force, in kJ/mol.");
    afir_energy_allowance.setDefaultValue(afirBase.energyAllowance);
    this->_fields.push_back(AFIROptimizerBase::afirEnergyAllowanceKey, afir_energy_allowance);

    UniversalSettings::IntDescriptor afir_phase_in(
        "The number of steps over which the full attractive force is slowly applied");
    afir_phase_in.setDefaultValue(afirBase.phaseIn);
    this->_fields.push_back(AFIROptimizerBase::afirPhaseInKey, afir_phase_in);

    UniversalSettings::BoolDescriptor afir_transfrom_coordinates(
        "Switch to transform the coordinates from Cartesian into an internal space.");
    afir_transfrom_coordinates.setDefaultValue(afirBase.transformCoordinates);
    this->_fields.push_back(AFIROptimizerBase::afirTransfromCoordinatesKey, afir_transfrom_coordinates);

    this->resetToDefaults();
  }
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_AFIROPTIMIZERSETTINGS_H_
