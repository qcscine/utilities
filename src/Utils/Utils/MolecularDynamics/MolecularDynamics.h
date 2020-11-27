/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MOLECULARDYNAMICS_H
#define UTILS_MOLECULARDYNAMICS_H

#include <Core/Interfaces/Calculator.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/MolecularTrajectory.h>
#include <Utils/Settings.h>

namespace Scine {

namespace Core {
struct Log;
}

namespace Utils {
class MDIntegrator;

class MolecularDynamics {
 public:
  /**
   * @brief Constructor that takes in a reference to a calculator.
   */
  explicit MolecularDynamics(Core::Calculator& calculator);
  /**
   * @brief Performs an MD simulation with the set calculator and the current settings.
   * @param structure The initial molecular structure.
   */
  void performMDSimulation(const AtomCollection& structure, Core::Log& log);
  /**
   * @brief Getter for the molecular trajectory.
   *
   *        The molecular trajectory holds the structures and energies along the MD simulation.
   *        The first element is always the initial structure with its energy.
   */
  Utils::MolecularTrajectory getMolecularTrajectory() const;
  /**
   * @brief Getter for the velocities corresponding to the structures in the molecular trajectory of the MD simulation.
   */
  std::vector<Utils::DisplacementCollection> getVelocities() const;
  /**
   * @brief Sets initial velocities for the MD simulation.
   */
  void setInitialVelocities(const Utils::DisplacementCollection& velocities);
  /**
   * @brief Accessor for the settings.
   * @return Utils::Settings& The settings.
   */
  Utils::Settings& settings();
  /**
   * @brief Constant accessor for the settings.
   * @return const Utils::Settings& The settings.
   */
  const Utils::Settings& settings() const;

 private:
  // The settings
  std::unique_ptr<Utils::Settings> settings_;
  // Apply the settings
  void applySettings();
  // The calculator
  Core::Calculator& calculator_;
  // Returns a unique pointer to a position updater.
  std::unique_ptr<MDIntegrator> getIntegrator() const;
  // Integration time step in femtoseconds
  double timeStep_;
  // Position updater algorithm
  std::string positionUpdaterAlgorithm_;
  // Initial velocities;
  Utils::DisplacementCollection initialVelocities_;
  // Option to couple to temperature bath
  bool temperatureBath_;
  // Target temperature if coupling to temperature bath is desired
  double temperature_;
  // The temperature bath relaxation time in units of the chosen time step.
  double relaxationTimeFactor_;
  // Number of MD steps
  int numberOfSteps_;
  // The molecular trajectory corresponding to the MD simulation (contains structures and energies)
  // The first element is always the initial structure.
  Utils::MolecularTrajectory structures_;
  // Decides whether velocities are recorded
  bool saveVelocities_;
  // A vector containing the velocities of each step
  std::vector<Utils::DisplacementCollection> velocities_;
  // The frequency with which the structures of the MD simulations and their energies
  // are saved to the molecular trajectory.
  int recordFrequency_;
  // Number of steps after which to remove the linear and angular momentum, respectively
  // If negative, do not remove.
  int linearMomentumRemovalFrequency_, angularMomentumRemovalFrequency_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_MOLECULARDYNAMICS_H
