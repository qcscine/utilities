/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "VelocityVerletMD.h"

namespace Scine {
namespace Utils {

Utils::DisplacementCollection VelocityVerletMD::calculateDisplacements(const Utils::GradientCollection& gradients) {
  assert(static_cast<int>(masses_.size()) == gradients.rows());
  previousAccelerations_ = accelerations_;
  calculateAccelerationsFromGradients(gradients);

  Utils::DisplacementCollection displacements = timeStep_ * (velocities_ + timeStep_ * accelerations_);
  velocities_ += (0.5 * timeStep_) * (accelerations_ + previousAccelerations_);
  rescaleVelocitiesForTemperatureBath();
  return displacements;
}

} // namespace Utils
} // namespace Scine