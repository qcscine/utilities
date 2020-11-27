/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_VELOCITYVERLETMD_H
#define UTILS_VELOCITYVERLETMD_H

#include "MDIntegrator.h"

namespace Scine {
namespace Utils {

/**
 * @class VelocityVerletMD VelocityVerletMD.h
 * @brief Class implementing the Velocity Verlet algorithm for Molecular Dynamics.
 *
 * The position is updated with
 * \f$ r_i(t + \Delta t) = r_i(t) + \Delta t v_i(t) + \frac{\Delta t^2}{2m_i} f_i(t) + O(\Delta t^3) \f$
 * and the velocities with
 * \f$ v_i(t + \Delta t) = v_i(t) + \frac{\Delta t}{2 m_i} (f_i(t) + f_i(t + \Delta t)) + O(\Delta t^3) \f$
 *
 */
class VelocityVerletMD : public MDIntegrator {
  /**
   * @brief Calculates the displacements from the gradients.
   */
  Utils::DisplacementCollection calculateDisplacements(const Utils::GradientCollection& gradients) override;
  /**
   * @brief The previous accelerations.
   */
  Utils::DisplacementCollection previousAccelerations_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_VELOCITYVERLETMD_H