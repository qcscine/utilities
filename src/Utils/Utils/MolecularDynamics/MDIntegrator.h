/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MDPOSITIONUPDATER_H
#define UTILS_MDPOSITIONUPDATER_H

#include <Utils/Typenames.h>

namespace Scine {
namespace Utils {

/**
 * @class MDPositionUpdater MDPositionUpdater.h
 * @brief This class manages the velocities and accelerations during an MD run and provides the displacements for
 *        a given set of gradients.
 */
class MDIntegrator {
 public:
  /// @brief Constructor that sets the default time step of 1 fs.
  MDIntegrator();
  /**
   * @brief Calculates the displacements from the gradients.
   */
  virtual Utils::DisplacementCollection calculateDisplacements(const Utils::GradientCollection& gradients) = 0;
  /**
   * @brief Sets the masses, initializes the velocities and accelerations to zero.
   * @param elements The element types of the molecular system.
   */
  void setElementTypes(const Utils::ElementTypeCollection& elements);
  /// @brief Resets velocities to zero.
  void resetVelocities();
  /**
   * @brief Sets the temperature around which the Maxwell-Boltzmann distribution for the initial velocities is centered.
   * @param T The temperature.
   */
  void setGenerationTemperatureInKelvin(double T);
  /// @brief Method to sample initial velocities from a Maxwell-Boltzmann distribution centered around the generation
  /// temperature
  void sampleVelocitiesFromBoltzmannDistribution();
  /// @brief Setter for the velocities.
  void setVelocities(const Utils::DisplacementCollection& velocities);
  /// @brief Getter for the velocities.
  Utils::DisplacementCollection getVelocities() const;
  /// @brief Calculates the current temperature.
  double getCurrentTemperature();
  /// @brief Sets the time step in the unit femtoseconds.
  void setTimeStepInFemtoseconds(double fs);
  /// @brief Sets the thermostat algorithm.
  void setThermostatAlgorithm(std::string thermostatAlgorithm);
  /// @brief Sets the thermostat time parameter in femtoseconds.
  void setTemperatureCouplingTimeInFemtoseconds(double fs);
  /// @brief Sets the target temperature in Kelvin.
  void setTargetTemperatureInKelvin(double T);
  //! @brief Sets the seed for noise generation in stochastic dynamics
  void setStochasticDynamicsSeed(int seed);
  /**
   * @brief Checks whether the integrator and thermostat algortihm are compatible.
   * @return true If the thermostat is available.
   * @return false If the thermostat is not available.
   */
  virtual bool checkThermostatAlgorithm() const = 0;
  //! @brief Sets a seed for the random initial velocities.
  void setSeed(int seed);
  //! @brief Removes the center of mass motion.
  void removeCenterOfMassLinearMomentum();
  //! @brief Removes the center of mass motion.
  void removeCenterOfMassAngularMomentum(const Eigen::MatrixX3d& positions);

  virtual ~MDIntegrator() = default;

 protected:
  //! @brief Updates the accelerations according to the given gradients.
  void calculateAccelerationsFromGradients(const Utils::GradientCollection& gradients);
  /**
   * @brief Apply Berendsen velocity scaling.
   */
  void rescaleVelocitiesWithBerendsen();
  std::vector<double> masses_;
  Utils::DisplacementCollection velocities_;
  Utils::DisplacementCollection accelerations_;
  int nAtoms_;      // number of atoms
  double timeStep_; // unit of this is not a.u., because masses are in u and not in terms of electron mass
  int seed_ = 42;   // seed for the initial velocities generation.

 protected:
  std::string thermostatAlgorithm_; // the thermostat algorithm
  double temperatureCouplingTime_;  // tau in a.u * sqrt(u / m_e)
  double targetTemperature_;        // k_B T in E_h
  int stochasticDynamicsSeed_;      // Seed used for noise generation in stochastic dynamics

 private:
  // Resets accelerations to zero.
  void resetAccelerations();
  // k_B T in E_h around which the initial Maxwell-Boltzmann distribution is centered
  double generationTemperature_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_MDPOSITIONUPDATER_H
