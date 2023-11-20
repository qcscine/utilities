/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MolecularDynamics.h"
#include "EulerMD.h"
#include "LeapFrogMD.h"
#include "MDIntegrator.h"
#include "MolecularDynamicsSettings.h"
#include "StochasticDynamics.h"
#include "VelocityVerletMD.h"
#include <Core/Log.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/GeometryUtilities.h>

namespace Scine {
namespace Utils {

MolecularDynamics::MolecularDynamics(Core::Calculator& calculator) : calculator_(calculator) {
  this->settings_ = std::make_unique<MolecularDynamicsSettings>();
  applySettings();
}

MolecularDynamics::MolecularDynamics(Core::CalculatorWithReference& calculatorWithReference)
  : calculatorWithReference_(calculatorWithReference), calculator_(calculatorWithReference.getReferenceCalculator()) {
  this->settings_ = std::make_unique<MolecularDynamicsSettings>();
  applySettings();
}

void MolecularDynamics::applySettings() {
  if (settings_->valid()) {
    generationSeed_ = settings_->getInt(SettingsNames::generationSeed);
    timeStep_ = settings_->getDouble(SettingsNames::timeStepInFemtoseconds);
    positionUpdaterAlgorithm_ = settings_->getString(SettingsNames::integrationAlgorithm);
    generationTemperature_ = settings_->getDouble(SettingsNames::generationTemperature);
    thermostatAlgorithm_ = settings_->getString(SettingsNames::thermostatAlgorithm);
    targetTemperature_ = settings_->getDouble(SettingsNames::targetTemperature);
    if (targetTemperature_ == 0.) {
      targetTemperature_ = generationTemperature_;
    }
    temperatureCouplingTime_ = settings_->getDouble(SettingsNames::temperatureCouplingTime);
    // Use correct temperature coupling defaults
    if (temperatureCouplingTime_ == 0.) {
      if (thermostatAlgorithm_ == OptionNames::berendsenThermostatOption) {
        temperatureCouplingTime_ = TemperatureCouplingDefaults::berendsenTimeDefault;
      }
      else if (positionUpdaterAlgorithm_ == OptionNames::stochasticDynamicsOption) {
        temperatureCouplingTime_ = TemperatureCouplingDefaults::stochasticDynamicsTimeDefault;
      }
    }
    stochasticDynamicsSeed_ = settings_->getInt(SettingsNames::stochasticDynamicsSeed);
    numberOfSteps_ = settings_->getInt(SettingsNames::numberOfMDSteps);
    recordFrequency_ = settings_->getInt(SettingsNames::recordFrequency);
    linearMomentumRemovalFrequency_ = settings_->getInt(SettingsNames::linearMomentumRemovalFrequency);
    angularMomentumRemovalFrequency_ = settings_->getInt(SettingsNames::angularMomentumRemovalFrequency);
    saveVelocities_ = settings_->getBool(SettingsNames::saveVelocities);
    saveTemperatures_ = settings_->getBool(SettingsNames::saveTemperatures);
  }
  else {
    settings_->throwIncorrectSettings();
  }
}

void MolecularDynamics::setCalculatorProperties() {
  Utils::PropertyList requiredProperties = Utils::Property::Energy | Utils::Property::Gradients;
  if (!calculator_.possibleProperties().containsSubSet(Utils::Property::Gradients)) {
    throw std::runtime_error("Gradients required for MD, but chosen calculator does not provide them.");
  }

  if (settings_->getBool(SettingsNames::requireCharges)) {
    if (!calculator_.possibleProperties().containsSubSet(Utils::Property::AtomicCharges)) {
      throw std::runtime_error("Atomic charges required, but chosen calculator does not provide them.");
    }
    requiredProperties.addProperty(Utils::Property::AtomicCharges);
  }

  if (settings_->getBool(SettingsNames::requireBondOrders)) {
    if (!calculator_.possibleProperties().containsSubSet(Utils::Property::BondOrderMatrix)) {
      throw std::runtime_error("Bond orders required, but chosen calculator does not provide them.");
    }
    requiredProperties.addProperty(Utils::Property::BondOrderMatrix);
  }
  calculator_.setRequiredProperties(requiredProperties);
}

Utils::MolecularTrajectory MolecularDynamics::getMolecularTrajectory() const {
  return structures_;
}

void MolecularDynamics::setInitialVelocities(const Utils::DisplacementCollection& velocities) {
  initialVelocities_ = velocities;
}

const Utils::Settings& MolecularDynamics::settings() const {
  return *settings_;
}

Utils::Settings& MolecularDynamics::settings() {
  return *settings_;
}

void MolecularDynamics::performMDSimulation(const AtomCollection& structure, Core::Log& log) {
  applySettings();
  setCalculatorProperties();
  structures_.clear();
  structures_.setElementTypes(structure.getElements());
  velocities_.clear();
  temperatures_.clear();

  // Initialize the position updater
  auto positionUpdater = getIntegrator();
  positionUpdater->setElementTypes(structure.getElements()); // Sets elements and resets velocities and accelerations
  positionUpdater->setSeed(generationSeed_);
  positionUpdater->setTimeStepInFemtoseconds(timeStep_);
  // Set target temperature and thermostat settings
  positionUpdater->setTargetTemperatureInKelvin(targetTemperature_);
  positionUpdater->setTemperatureCouplingTimeInFemtoseconds(temperatureCouplingTime_);
  positionUpdater->setThermostatAlgorithm(thermostatAlgorithm_);
  if (!positionUpdater->checkThermostatAlgorithm()) {
    throw Core::InitializationException("The chosen thermostat algorithm and MD integrator are incompatible!");
  }
  if (positionUpdaterAlgorithm_ == OptionNames::stochasticDynamicsOption) {
    positionUpdater->setStochasticDynamicsSeed(stochasticDynamicsSeed_);
  }

  // Set or generate initial velocities
  if (initialVelocities_.size() != 0) {
    positionUpdater->setVelocities(initialVelocities_);
  }
  else if (generationTemperature_ > 0.) {
    positionUpdater->setGenerationTemperatureInKelvin(generationTemperature_);
    positionUpdater->sampleVelocitiesFromBoltzmannDistribution();
  }

  // Initialize the calculator
  calculator_.setStructure(structure);

  // First calculation
  const Results& initialResults = this->calculateWithCorrectCalculator();

  if (!initialResults.get<Property::SuccessfulCalculation>()) {
    throw std::runtime_error("Initial gradient calculation in MD simulation failed.");
  }
  double energy = initialResults.get<Property::Energy>();
  auto gradients = initialResults.get<Property::Gradients>();
  Utils::PositionCollection positions = structure.getPositions();

  // Apply bias potential
  if (biasPotential_) {
    auto bias = biasPotential_(positions, initialResults, 0);
    gradients += bias.second;
  }

  structures_.push_back(positions, energy);
  Utils::DisplacementCollection displacements;
  if (saveVelocities_) {
    velocities_.push_back(positionUpdater->getVelocities());
  }
  if (saveTemperatures_) {
    temperatures_.push_back(positionUpdater->getCurrentTemperature() * Constants::joule_per_hartree / Constants::boltzmannConstant);
  }

  log.output << "Calculation for initial structure done." << Core::Log::nl;
  log.output << "Starting MD simulation with " << numberOfSteps_ << " steps..." << Core::Log::nl;

  // MD simulation
  for (int i = 0; i < numberOfSteps_; ++i) {
    displacements = positionUpdater->calculateDisplacements(gradients);
    // Update Positions
    positions += displacements;
    // Remove center of mass linear and angular momentum
    if (i > 0 && (i % linearMomentumRemovalFrequency_) == 0) {
      positionUpdater->removeCenterOfMassLinearMomentum();
    }
    if (i > 0 && (i % angularMomentumRemovalFrequency_) == 0) {
      positionUpdater->removeCenterOfMassAngularMomentum(positions);
    }

    calculator_.modifyPositions(positions);

    const Results& results = this->calculateWithCorrectCalculator();

    if (!results.get<Property::SuccessfulCalculation>()) {
      finalVelocities_ = positionUpdater->getVelocities();
      throw std::runtime_error("Gradient calculation in MD simulation failed.");
    }

    gradients = results.get<Property::Gradients>();

    // Apply bias potential
    if (biasPotential_) {
      auto bias = biasPotential_(positions, results, i + 1);
      gradients += bias.second;
    }

    if ((i + 1) % recordFrequency_ == 0) {
      structures_.push_back(positions, results.get<Property::Energy>());
      if (saveVelocities_) {
        velocities_.push_back(positionUpdater->getVelocities());
      }
      if (saveTemperatures_) {
        // Transform from kbT in Hartree to T in Kelvin and save
        temperatures_.push_back(positionUpdater->getCurrentTemperature() * Constants::joule_per_hartree /
                                Constants::boltzmannConstant);
      }
      log.debug << "The structure number " << ((i + 1) / recordFrequency_) + 1
                << " has been added to the MD trajectory." << Core::Log::nl;
    }
    if (externalStop_) {
      auto shouldStop = externalStop_(positions, results, i + 1);
      if (shouldStop) {
        log.output.line("External stop function demands termination of MD.");
        break;
      }
    }
  }
  finalVelocities_ = positionUpdater->getVelocities();
  log.output.line("MD simulation done.");
}

const Utils::Results& MolecularDynamics::calculateWithCorrectCalculator() const {
  if (calculatorWithReference_) {
    return calculatorWithReference_->calculate();
  }
  return calculator_.calculate("");
}

std::unique_ptr<MDIntegrator> MolecularDynamics::getIntegrator() const {
  if (positionUpdaterAlgorithm_ == OptionNames::leapFrogOption) {
    return std::make_unique<LeapFrogMD>();
  }

  if (positionUpdaterAlgorithm_ == OptionNames::eulerOption) {
    return std::make_unique<EulerMD>();
  }
  if (positionUpdaterAlgorithm_ == OptionNames::stochasticDynamicsOption) {
    return std::make_unique<StochasticDynamics>();
  }
  assert(positionUpdaterAlgorithm_ == OptionNames::velocityVerletOption);

  return std::make_unique<VelocityVerletMD>();
}

std::vector<Utils::DisplacementCollection> MolecularDynamics::getVelocities() const {
  return velocities_;
}

Utils::DisplacementCollection MolecularDynamics::getFinalVelocities() const {
  return finalVelocities_;
}

std::vector<double> MolecularDynamics::getTemperatures() const {
  return temperatures_;
}

void MolecularDynamics::setExternalStop(std::function<bool(const PositionCollection&, const Results&, int)> externalStop) {
  externalStop_ = std::move(externalStop);
}

void MolecularDynamics::setBiasPotential(
    std::function<std::pair<double, GradientCollection>(const PositionCollection&, const Results&, int)> biasPotential) {
  biasPotential_ = std::move(biasPotential);
}

} // namespace Utils
} // namespace Scine
