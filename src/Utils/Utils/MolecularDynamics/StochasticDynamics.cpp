/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "StochasticDynamics.h"
#include "MolecularDynamicsSettings.h"
#include <random>

namespace Scine {
namespace Utils {

bool StochasticDynamics::checkThermostatAlgorithm() const {
  return thermostatAlgorithm_ == OptionNames::noThermostatOption;
}

Utils::DisplacementCollection StochasticDynamics::calculateDisplacements(const Utils::GradientCollection& gradients) {
  if (!preparedScaling_) {
    prepareScaling();
    preparedScaling_ = true;
  }
  assert(static_cast<int>(masses_.size()) == gradients.rows());
  calculateAccelerationsFromGradients(gradients);

  // Usual leapfrog velocity update (1)
  velocities_ += timeStep_ * accelerations_;
  // Friction application (2)
  createNoise();
  // noise and noiseScaling are arrays i.e. multiplied elementwise
  Utils::DisplacementCollection scaledVelocities = -velocityScaling_ * velocities_ + (noiseScaling_ * noise_).matrix();
  // Displacements: half step with usual v and half step with friction scaling (3)
  Utils::DisplacementCollection displacements = timeStep_ * (velocities_ + 0.5 * scaledVelocities);
  // Velocity update (4)
  velocities_ += scaledVelocities;

  return displacements;
}

void StochasticDynamics::prepareScaling() {
  // Seed the random number generator for noise production
  gen_.seed(stochasticDynamicsSeed_);
  // Compute velocity prefactor assuming that the time coupling parameter is the inverse of the friction constant
  velocityScaling_ = 1 - std::exp(-timeStep_ / temperatureCouplingTime_);
  // Compute scaling array for noise
  noiseScaling_.setConstant(nAtoms_, 3, velocityScaling_ * (2 - velocityScaling_) * targetTemperature_);
  for (int i = 0; i != nAtoms_; ++i) {
    noiseScaling_.row(i) /= masses_[i];
  }
  noiseScaling_ = noiseScaling_.sqrt();
}

void StochasticDynamics::createNoise() {
  std::normal_distribution<> d(0., 1.);
  noise_ = Eigen::ArrayXXd::NullaryExpr(nAtoms_, 3, [&]() { return d(gen_); });
}

} // namespace Utils
} // namespace Scine
