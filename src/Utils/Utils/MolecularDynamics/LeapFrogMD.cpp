/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LeapFrogMD.h"
#include "MolecularDynamicsSettings.h"

namespace Scine {
namespace Utils {

bool LeapFrogMD::checkThermostatAlgorithm() const {
  return thermostatAlgorithm_ == OptionNames::berendsenThermostatOption ||
         thermostatAlgorithm_ == OptionNames::noThermostatOption;
}

Utils::DisplacementCollection LeapFrogMD::calculateDisplacements(const Utils::GradientCollection& gradients) {
  assert(static_cast<int>(masses_.size()) == gradients.rows());
  calculateAccelerationsFromGradients(gradients);

  // Velocities one half step ahead
  velocities_ += timeStep_ * accelerations_;
  if (thermostatAlgorithm_ == OptionNames::berendsenThermostatOption) {
    rescaleVelocitiesWithBerendsen();
  }
  Utils::DisplacementCollection displacements = timeStep_ * velocities_;
  return displacements;
}

} // namespace Utils
} // namespace Scine
