/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_LEAPFROGMD_H
#define UTILS_LEAPFROGMD_H

#include "MDIntegrator.h"

namespace Scine {
namespace Utils {

/**
 * @class LeapFrogMD LeapFrogMD.h
 * @brief Class implementing the LeapFrog algorithm for Molecular Dynamics.
 *
 * The velocities are updated with
 * \f$ v_i(t + \frac{1}{2} \Delta t) = v_i(t - \frac{1}{2} \Delta t) + \frac{\Delta t}{m_i} f_i(t) \f$
 * and the positions with
 * \f$ r_i(t + \Delta t) = r_i(t) + \Delta t v_i(t + \frac{1}{2} \Delta t) \f$
 *
 */
class LeapFrogMD : public MDIntegrator {
  /**
   * @brief Calculates the displacements from the gradients.
   */
  Utils::DisplacementCollection calculateDisplacements(const Utils::GradientCollection& gradients) override;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_LEAPFROGMD_H
