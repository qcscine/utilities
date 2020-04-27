/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LeapFrogMD.h"

namespace Scine {
namespace Utils {

Utils::DisplacementCollection LeapFrogMD::calculateDisplacements(const Utils::GradientCollection& gradients) {
  assert(static_cast<int>(masses_.size()) == gradients.rows());
  calculateAccelerationsFromGradients(gradients);

  // Velocities one half step ahead
  velocities_ += timeStep_ * accelerations_;
  rescaleVelocitiesForTemperatureBath();
  Utils::DisplacementCollection displacements = timeStep_ * velocities_;
  return displacements;
}

} // namespace Utils
} // namespace Scine