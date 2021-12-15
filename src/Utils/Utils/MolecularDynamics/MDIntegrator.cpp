/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MDIntegrator.h"
#include <Utils/Constants.h>
#include <Utils/Geometry/GeometryUtilities.h>
#include <Eigen/Geometry>
#include <random>

namespace Scine {
namespace Utils {

namespace {
constexpr double femtosecondsToAtomicUnits =
    2. * Constants::pi * Constants::joule_per_hartree * 1.e-15 / Constants::planckConstant; // 1 fs in atomic time units
constexpr double kelvinPerHartree = Constants::joule_per_hartree / Constants::boltzmannConstant; // k_B T in E_h to T in K
constexpr double hartreePerKelvin = 1. / kelvinPerHartree; // T in K to k_B T in E_h
} // namespace

MDIntegrator::MDIntegrator()
  : targetTemperature_(300. * hartreePerKelvin), generationTemperature_(300. * hartreePerKelvin) {
  setTimeStepInFemtoseconds(1.0);
  setTemperatureCouplingTimeInFemtoseconds(10.0);
}

void MDIntegrator::setElementTypes(const Utils::ElementTypeCollection& elements) {
  nAtoms_ = static_cast<int>(elements.size());
  masses_ = Utils::Geometry::Properties::getMasses(elements);
  resetVelocities();
  resetAccelerations();
}

void MDIntegrator::calculateAccelerationsFromGradients(const Utils::GradientCollection& gradients) {
  assert(gradients.size() == accelerations_.size());
  for (int i = 0; i < nAtoms_; ++i) {
    accelerations_.row(i) = (-1.0 / masses_[i]) * gradients.row(i);
  }
}

void MDIntegrator::rescaleVelocitiesWithBerendsen() {
  double currentTemperature = getCurrentTemperature();
  double factor = std::sqrt(1. + timeStep_ / temperatureCouplingTime_ * (targetTemperature_ / currentTemperature - 1));
  velocities_ *= factor;
}

void MDIntegrator::resetVelocities() {
  velocities_.resize(nAtoms_, 3);
  velocities_.setZero();
}

void MDIntegrator::resetAccelerations() {
  accelerations_.resize(nAtoms_, 3);
  accelerations_.setZero();
}

void MDIntegrator::setVelocities(const Utils::DisplacementCollection& velocities) {
  velocities_ = velocities;
}

void MDIntegrator::setGenerationTemperatureInKelvin(double T) {
  generationTemperature_ = T * hartreePerKelvin;
}

void MDIntegrator::sampleVelocitiesFromBoltzmannDistribution() {
  std::mt19937 gen(seed_);
  double numerator = std::sqrt(generationTemperature_);
  int index = 0;
  for (double mass : masses_) {
    std::normal_distribution<> d(0., numerator * std::sqrt(1. / mass));
    velocities_.row(index) = Position{d(gen), d(gen), d(gen)};
    ++index;
  }
}

Utils::DisplacementCollection MDIntegrator::getVelocities() const {
  return velocities_;
}

void MDIntegrator::setTimeStepInFemtoseconds(double fs) {
  timeStep_ = fs * femtosecondsToAtomicUnits / sqrt(Constants::electronRestMass_per_u);
}

void MDIntegrator::setThermostatAlgorithm(std::string thermostatAlgorithm) {
  thermostatAlgorithm_ = thermostatAlgorithm;
}

void MDIntegrator::setTemperatureCouplingTimeInFemtoseconds(double fs) {
  temperatureCouplingTime_ = fs * femtosecondsToAtomicUnits / sqrt(Constants::electronRestMass_per_u);
}

void MDIntegrator::setTargetTemperatureInKelvin(double T) {
  targetTemperature_ = T * hartreePerKelvin;
}

double MDIntegrator::getCurrentTemperature() {
  int numberDegreesOfFreedom = 3 * nAtoms_; // TODO: Subtract translation and remove center of mass motion
  double currentTemperature =
      (velocities_.rowwise().squaredNorm().array() * Eigen::Map<const Eigen::ArrayXd>(masses_.data(), masses_.size())).sum() /
      numberDegreesOfFreedom;
  return currentTemperature;
}

void MDIntegrator::setSeed(int seed) {
  seed_ = seed;
}

void MDIntegrator::setStochasticDynamicsSeed(int seed) {
  stochasticDynamicsSeed_ = seed;
}

void MDIntegrator::removeCenterOfMassLinearMomentum() {
  const Eigen::VectorXd& masses = Eigen::Map<const Eigen::VectorXd>(masses_.data(), masses_.size());
  double totalMass = masses.sum();

  // Remove total linear momentum
  Eigen::MatrixX3d linearMomentumVector = velocities_.array().colwise() * masses.array();
  Eigen::RowVector3d centerOfMassVelocity = linearMomentumVector.colwise().sum() / totalMass;

  velocities_.rowwise() -= centerOfMassVelocity;
}

void MDIntegrator::removeCenterOfMassAngularMomentum(const Eigen::MatrixX3d& positions) {
  const Eigen::VectorXd& masses = Eigen::Map<const Eigen::VectorXd>(masses_.data(), masses_.size());
  auto centerOfMass = Utils::Geometry::Properties::getCenterOfMass(positions, masses_);
  auto inertiaTensor = Utils::Geometry::Properties::calculateInertiaTensor(positions, masses_, centerOfMass);
  Eigen::MatrixX3d positionsRelativeToCOM = positions.rowwise() - centerOfMass;

  // Remove total angular momentum
  Eigen::MatrixX3d linearMomentumVector = velocities_.array().colwise() * masses.array();
  Eigen::RowVector3d totalAngularMomentum = Eigen::RowVector3d::Zero(1, 3);
  for (int i = 0; i < linearMomentumVector.rows(); ++i) {
    Eigen::RowVector3d particleAngularMomentum = positionsRelativeToCOM.row(i).cross(linearMomentumVector.row(i));
    totalAngularMomentum += particleAngularMomentum;
  }

  Eigen::RowVector3d angularVelocityVector =
      inertiaTensor.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(totalAngularMomentum.transpose());

  for (int i = 0; i < positionsRelativeToCOM.rows(); ++i) {
    velocities_.row(i) -= angularVelocityVector.cross(positionsRelativeToCOM.row(i));
  }
}
} // namespace Utils
} // namespace Scine
