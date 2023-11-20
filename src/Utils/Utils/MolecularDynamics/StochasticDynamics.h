/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_STOCHASTICDYNAMICS_H
#define UTILS_STOCHASTICDYNAMICS_H

#include "MDIntegrator.h"
#include <random>

namespace Scine {
namespace Utils {

/**
 * @class StochasticDynamics StochasticDynamics.h
 * @brief Class implementing stochastic dynamics.
 *
 * A leapfrog-based stochastic dynamics integrator:
 *
 * \f$ v' = v(t- \frac{1}{2}\Delta t) + a \Delta t \f$ <BR>
 * \f$ \Delta v = - f v' + \sqrt{f(2-f)(\frac{k_B T_{target}}{m})}\xi\f$ <BR>
 * \f$ x(t + \Delta t) = x(t) + (v' + \frac{1}{2}\Delta v) \Delta t\f$ <BR>
 * \f$ v(t + \frac{1}{2} \Delta t) = v' + \Delta v\f$ <BR>
 * \f$ f = 1 - \exp(-\gamma \Delta t) = 1 - \exp(-\frac{\Delta t}{\tau})\f$ <BR>
 *
 * &xi; is drawn from a normal distribution. <BR>
 * This integrator always employs the target temperature. It cannot be coupled to a separate thermostat.
 * The friction constant &gamma; is set to the inverse of the temperature coupling time parameter.
 *
 *
 * Implemented after
 * J. Chem. Theory Comput. 2012, 8 (10), 3637–3649.
 * https://doi.org/10.1021/ct3000876
 *
 */
class StochasticDynamics : public MDIntegrator {
 public:
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

 private:
  // Whether the scaling prefactors were initialized
  bool preparedScaling_ = false;
  // Random number generator for noise generation
  std::mt19937 gen_;
  // Velocity scaling factor encorporating friction
  double velocityScaling_;
  // Noise scaling encorporating friction and target temperature
  Eigen::ArrayXXd noiseScaling_;
  // Gaussian distributed noise array
  Eigen::ArrayXXd noise_;
  /// @brief Computes the scaling factors
  void prepareScaling();
  /// @brief Updates the noise array
  void createNoise();
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_STOCHASTICDYNAMICS_H
