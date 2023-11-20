/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "EulerMD.h"
#include "MolecularDynamicsSettings.h"

namespace Scine {
namespace Utils {

bool EulerMD::checkThermostatAlgorithm() const {
  return thermostatAlgorithm_ == OptionNames::berendsenThermostatOption ||
         thermostatAlgorithm_ == OptionNames::noThermostatOption;
}

Utils::DisplacementCollection EulerMD::calculateDisplacements(const Utils::GradientCollection& gradients) {
  assert(static_cast<int>(masses_.size()) == gradients.rows());
  calculateAccelerationsFromGradients(gradients);

  Utils::DisplacementCollection displacements = timeStep_ * (velocities_ + (0.5 * timeStep_) * accelerations_);
  velocities_ += timeStep_ * accelerations_;
  if (thermostatAlgorithm_ == OptionNames::berendsenThermostatOption) {
    rescaleVelocitiesWithBerendsen();
  }
  return displacements;
}

} // namespace Utils
} // namespace Scine
