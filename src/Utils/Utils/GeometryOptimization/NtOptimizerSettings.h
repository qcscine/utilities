/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
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
    this->_fields.push_back(NtOptimizer::ntSdFactorKey, factor);

    UniversalSettings::DoubleDescriptor repulsiveStop(
        "The stop parameter given in multiples/fractions of covalent radii sums in the repulsive case.");
    repulsiveStop.setDefaultValue(nt.check.repulsiveStop);
    this->_fields.push_back(NtOptimizer::ntRepulsiveStopKey, repulsiveStop);

    UniversalSettings::DoubleDescriptor attractiveStop(
        "The stop parameter given in multiples/fractions of covalent radii sums in the attractive case.");
    attractiveStop.setDefaultValue(nt.check.attractiveStop);
    this->_fields.push_back(NtOptimizer::ntAttractiveStopKey, attractiveStop);

    UniversalSettings::IntDescriptor max_iter("The maximum number of iterations.");
    max_iter.setDefaultValue(nt.check.maxIter);
    this->_fields.push_back(NtOptimizer::ntMaxIterKey, max_iter);

    UniversalSettings::IntListDescriptor nt_rhs_list(
        "The list indices of atoms to be artificially forced onto or away from those in the LHS list.");
    nt_rhs_list.setDefaultValue(nt.rhsList);
    this->_fields.push_back(NtOptimizer::ntRHSListKey, nt_rhs_list);

    UniversalSettings::IntListDescriptor nt_lhs_list(
        "The list indices of atoms to be artificially forced onto or away from those in the RHS list.");
    nt_lhs_list.setDefaultValue(nt.lhsList);
    this->_fields.push_back(NtOptimizer::ntLHSListKey, nt_lhs_list);

    UniversalSettings::BoolDescriptor nt_attractive("Switch for the artificial force to be attractive or repulsive.");
    nt_attractive.setDefaultValue(nt.attractive);
    this->_fields.push_back(NtOptimizer::ntAttractiveKey, nt_attractive);

    UniversalSettings::DoubleDescriptor nt_total_force_norm(
        "The norm of the summed additional forces acting on all listed atoms.");
    nt_total_force_norm.setDefaultValue(nt.totalForceNorm);
    this->_fields.push_back(NtOptimizer::ntTotalForceNormKey, nt_total_force_norm);

    UniversalSettings::BoolDescriptor nt_transform_coordinates(
        "Switch to transform the coordinates from Cartesian into an internal space.");
    nt_transform_coordinates.setDefaultValue(nt.transformCoordinates);
    this->_fields.push_back(NtOptimizer::ntTransfromCoordinatesKey, nt_transform_coordinates);

    this->resetToDefaults();
  }
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_NTOPTIMIZERSETTINGS_H_
