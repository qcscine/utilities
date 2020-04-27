/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MolecularDynamics.h"
#include "EulerMD.h"
#include "LeapFrogMD.h"
#include "MDIntegrator.h"
#include "MolecularDynamicsSettings.h"
#include "VelocityVerletMD.h"
#include <Utils/CalculatorBasics.h>
#include <Utils/Geometry/GeometryUtilities.h>
#include <Utils/IO/Logger.h>

namespace Scine {
namespace Utils {

MolecularDynamics::MolecularDynamics(Core::Calculator& calculator) : calculator_(calculator) {
  calculator_.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  this->settings_ = std::make_unique<MolecularDynamicsSettings>();
  applySettings();
}

void MolecularDynamics::applySettings() {
  if (settings_->check()) {
    timeStep_ = settings_->getDouble(SettingsNames::timeStepInFemtoseconds);
    positionUpdaterAlgorithm_ = settings_->getString(SettingsNames::integrationAlgorithm);
    temperatureBath_ = settings_->getBool(SettingsNames::temperatureBath);
    temperature_ = settings_->getDouble(SettingsNames::targetTemperature);
    relaxationTimeFactor_ = settings_->getDouble(SettingsNames::relaxationTimeFactor);
    numberOfSteps_ = settings_->getInt(SettingsNames::numberOfMDSteps);
    recordFrequency_ = settings_->getInt(SettingsNames::recordFrequency);
    linearMomentumRemovalFrequency_ = settings_->getInt(SettingsNames::linearMomentumRemovalFrequency);
    angularMomentumRemovalFrequency_ = settings_->getInt(SettingsNames::angularMomentumRemovalFrequency);
    saveVelocities_ = settings_->getBool(SettingsNames::saveVelocities);
  }
  else {
    throw Core::InitializationException("Settings are invalid!");
  }
}

Utils::MolecularTrajectory MolecularDynamics::getMolecularTrajectory() const {
  return structures_;
};

void MolecularDynamics::setInitialVelocities(const Utils::DisplacementCollection& velocities) {
  initialVelocities_ = velocities;
}

const Utils::Settings& MolecularDynamics::settings() const {
  return *settings_;
}

Utils::Settings& MolecularDynamics::settings() {
  return *settings_;
}

void MolecularDynamics::performMDSimulation(const AtomCollection& structure) {
  applySettings();
  structures_.clear();
  structures_.setElementTypes(structure.getElements());
  velocities_.clear();

  // Initialize the position updater
  auto positionUpdater = getIntegrator();
  positionUpdater->setElementTypes(structure.getElements());
  positionUpdater->setTimeStepInFemtoseconds(timeStep_);
  positionUpdater->setRelaxationTimeFactor(relaxationTimeFactor_);
  if (initialVelocities_.size() != 0)
    positionUpdater->setVelocities(initialVelocities_);
  if (temperatureBath_)
    positionUpdater->setTargetTemperatureInKelvin(temperature_);

  // Initialize the calculator
  calculator_.setStructure(structure);

  // First calculation
  const auto& initialResults = calculator_.calculate("");
  double energy = initialResults.get<Property::Energy>();
  auto gradients = initialResults.get<Property::Gradients>();

  Utils::PositionCollection positions = structure.getPositions();
  Utils::DisplacementCollection displacements;
  structures_.push_back(positions, energy);
  if (saveVelocities_)
    velocities_.push_back(initialVelocities_);

  Utils::Log::info() << "Calculation for initial structure done.";
  Utils::Log::info() << "Starting MD simulation with " << numberOfSteps_ << " steps...";

  // MD simulation
  for (int i = 0; i < numberOfSteps_; ++i) {
    displacements = positionUpdater->calculateDisplacements(gradients);
    // Update Positions
    positions += displacements;
    // Remove center of mass linear and angular momentum
    if (i > 0 && (i % linearMomentumRemovalFrequency_) == 0)
      positionUpdater->removeCenterOfMassLinearMomentum(positions);
    if (i > 0 && (i % angularMomentumRemovalFrequency_) == 0)
      positionUpdater->removeCenterOfMassAngularMomentum(positions);

    calculator_.modifyPositions(positions);
    const auto& results = calculator_.calculate("");

    gradients = results.get<Property::Gradients>();
    if ((i + 1) % recordFrequency_ == 0) {
      structures_.push_back(positions, results.get<Property::Energy>());
      if (saveVelocities_)
        velocities_.push_back(positionUpdater->getVelocities());
      Utils::Log::debug() << "The structure number " << ((i + 1) / recordFrequency_) + 1
                          << " has been added to the MD trajectory.";
    }
  }
  Utils::Log::info() << "MD simulation done.";
}

std::unique_ptr<MDIntegrator> MolecularDynamics::getIntegrator() const {
  if (positionUpdaterAlgorithm_ == OptionNames::leapFrogOption)
    return std::make_unique<LeapFrogMD>();

  if (positionUpdaterAlgorithm_ == OptionNames::eulerOption)
    return std::make_unique<EulerMD>();

  assert(positionUpdaterAlgorithm_ == OptionNames::velocityVerletOption);

  return std::make_unique<VelocityVerletMD>();
}

std::vector<Utils::DisplacementCollection> MolecularDynamics::getVelocities() const {
  return velocities_;
}

} // namespace Utils
} // namespace Scine
