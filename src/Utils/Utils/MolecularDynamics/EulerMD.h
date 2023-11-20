/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_EULERMD_H
#define UTILS_EULERMD_H

#include "MDIntegrator.h"

namespace Scine {
namespace Utils {

/**
 * @class EulerMD EulerMD.h
 * @brief Class implementing the Euler algorithm for Molecular Dynamics.
 *
 * The position is updated with
 * \f$ r_i(t + \Delta t) = r_i(t) + \Delta t v_i(t) + \frac{\Delta t^2}{2m_i} f_i(t) + O(\Delta t^3) \f$
 * and the velocities with
 * \f$ v_i(t + \Delta t) = v_i(t) + \frac{\Delta t}{m_i} f_i(t) + O(\Delta t^2) \f$
 *
 * This integrator can be combined with the Berendsen thermostat.
 */
class EulerMD : public MDIntegrator {
  /**
   * @brief Checks whether the set integrator and thermostat algortihm are compatible.
   * @return true If the thermostat is available.
   * @return false If the thermostat is not available.
   */
  bool checkThermostatAlgorithm() const override;

  /**
   * @brief Calculates the displacements from the gradients.
   */
  Utils::DisplacementCollection calculateDisplacements(const Utils::GradientCollection& gradients) override;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_EULERMD_H
