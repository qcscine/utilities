/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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
  /// @brief Setter for the velocities.
  void setVelocities(const Utils::DisplacementCollection& velocities);
  /// @brief Getter for the velocities.
  Utils::DisplacementCollection getVelocities() const;
  /// @brief Sets the time step in the unit femtoseconds.
  void setTimeStepInFemtoseconds(double fs);
  /// @brief The temperature bath relaxation time in units of the chosen time step.
  void setRelaxationTimeFactor(double factor);
  /// @brief Sets the target temperature in Kelvin.
  void setTargetTemperatureInKelvin(double T);
  //! @brief Sets a seed for the random initial velocities.
  void setSeed(int seed);
  //! @brief Removes the center of mass motion.
  void removeCenterOfMassLinearMomentum(const Eigen::MatrixX3d& positions);
  //! @brief Removes the center of mass motion.
  void removeCenterOfMassAngularMomentum(const Eigen::MatrixX3d& positions);

  virtual ~MDIntegrator() = default;

 protected:
  //! @brief Updates the accelerations according to the given gradients.
  void calculateAccelerationsFromGradients(const Utils::GradientCollection& gradients);
  //! @brief Berendsen thermostat.
  void rescaleVelocitiesForTemperatureBath();
  //! @brief Method to sample initial velocities from a Maxwell-Boltzmann distribution.
  void sampleVelocitiesFromBoltzmannDistribution();
  std::vector<double> masses_;
  Utils::DisplacementCollection velocities_;
  Utils::DisplacementCollection accelerations_;
  int nAtoms_;      // number of atoms
  double timeStep_; // unit of this is not a.u., because masses are in u and not in terms of electron mass
  int seed_ = 42;   // seed for the initial velocities generation.

 private:
  // Resets accelerations to zero.
  void resetAccelerations();
  // Calculates the current temperature.
  double getCurrentTemperature();
  // Option to use temperature bath is switched on when a target temperature is set.
  bool temperatureBath_ = false;
  double temperatureRelaxationTime_;
  // Default to room temperature for initial velocities.
  double targetTemperature_;
  double relaxationTimeFactor_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_MDPOSITIONUPDATER_H
