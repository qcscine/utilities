/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MOLECULARDYNAMICS_H
#define UTILS_MOLECULARDYNAMICS_H

#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/CalculatorWithReference.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/MolecularTrajectory.h>
#include <Utils/Typenames.h>
#include <boost/optional.hpp>

namespace Scine {

namespace Core {
struct Log;
}

namespace Utils {
class MDIntegrator;
class Results;
class Settings;
class AtomCollection;

class MolecularDynamics {
 public:
  /**
   * @brief Constructor that takes in a reference to a calculator.
   */
  explicit MolecularDynamics(Core::Calculator& calculator);
  /**
   * @brief Constructor that takes in a reference to a calculator with a reference calculator.
   */
  explicit MolecularDynamics(Core::CalculatorWithReference& calculatorWithReference);
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
   * @return std::vector<Utils::DisplacementCollection> The velocities.
   */
  std::vector<Utils::DisplacementCollection> getVelocities() const;
  /**
   * @brief Getter for the last velocities encountered during the MD simulation.
   * @note These velocities can be different from the last element of the velocities obtainable via getVelocities() if
   * the record frequency is not one.
   * @return Utils::DisplacementCollection The velocities.
   */
  Utils::DisplacementCollection getFinalVelocities() const;
  /**
   * @brief Getter for the temperatures corresponding to the structures in the molecular trajectory of the MD simulation.
   * @return const std::vector<double> The temperatures in Kelvin.
   */
  std::vector<double> getTemperatures() const;
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
  /**
   * @brief Setter for the external stop function.
   * @param Function with arguments: current positions, calculation results, number of steps since start of MD.
   */
  void setExternalStop(std::function<bool(const PositionCollection&, const Results&, int)> externalStop);
  /**
   * @brief Setter for the bias potential function.
   * @param Function with arguments: current positions, calculation results, number of steps since start of MD.
   */
  void setBiasPotential(std::function<std::pair<double, GradientCollection>(const PositionCollection&, const Results&, int)> biasPotential);

 private:
  // The settings
  std::unique_ptr<Utils::Settings> settings_;
  // Apply the settings
  void applySettings();
  // Sets the required calculator properties and checks whether all requested ones are available
  void setCalculatorProperties();
  // An optional calculator with reference calculator
  boost::optional<Core::CalculatorWithReference&> calculatorWithReference_;
  // The (reference) calculator
  Core::Calculator& calculator_;
  // Runs a calculation either with the calculator with reference if it exists or directly with the calculator if not
  const Utils::Results& calculateWithCorrectCalculator() const;
  // Returns a unique pointer to a position updater.
  std::unique_ptr<MDIntegrator> getIntegrator() const;
  // The seed for the initial velocity generation
  int generationSeed_;
  // Integration time step in femtoseconds
  double timeStep_;
  // Position updater algorithm
  std::string positionUpdaterAlgorithm_;
  // Initial velocities;
  Utils::DisplacementCollection initialVelocities_;
  // Thermostat algorithm
  std::string thermostatAlgorithm_;
  // Target temperature if coupling to temperature bath is desired
  double targetTemperature_;
  // Temperature for which initial velocities are drawn from a Boltzmann distribution if they are not given explicitly
  double generationTemperature_;
  // The thermostat coupling time in femtoseconds
  double temperatureCouplingTime_;
  // The seed used for stochastic dynamics
  int stochasticDynamicsSeed_;
  // Number of MD steps
  int numberOfSteps_;
  // The molecular trajectory corresponding to the MD simulation (contains structures and energies)
  // The first element is always the initial structure.
  Utils::MolecularTrajectory structures_;
  // Decides whether velocities are recorded
  bool saveVelocities_;
  // Decides whether temperatures are recorded
  bool saveTemperatures_;

  // A vector containing the velocities of each step
  std::vector<Utils::DisplacementCollection> velocities_;
  // The final velocities encountered duting the MD simulation
  // These can be different from velocities_.back() if the record frequency is not one
  Utils::DisplacementCollection finalVelocities_;
  // A vector containing the temperatures of each step in Kelvin
  std::vector<double> temperatures_;
  // The frequency with which the structures of the MD simulations and their energies
  // are saved to the molecular trajectory.
  int recordFrequency_;
  // Number of steps after which to remove the linear and angular momentum, respectively
  // If negative, do not remove.
  int linearMomentumRemovalFrequency_, angularMomentumRemovalFrequency_;
  // The external stop function.
  std::function<bool(const PositionCollection&, const Results&, int)> externalStop_ = nullptr;
  // The bias potential.
  std::function<std::pair<double, GradientCollection>(const PositionCollection&, const Results&, int)> biasPotential_ = nullptr;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_MOLECULARDYNAMICS_H
