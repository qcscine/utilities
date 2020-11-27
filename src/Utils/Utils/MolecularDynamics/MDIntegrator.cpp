/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MDIntegrator.h"
#include <Utils/Geometry/GeometryUtilities.h>
#include <Eigen/Geometry>
#include <random>

namespace Scine {
namespace Utils {

static constexpr double femtosecondsToAtomicUnits = 41.3413746; // 1 fs in atomic time units
static constexpr double massConversionFactor = 1822.88848619;   // 1 u in electron masses
static constexpr double atomicUnitsToKelvin = 3.1577464e5;
static constexpr double kelvinToAtomicUnits = 1.0 / atomicUnitsToKelvin;
static constexpr double boltzmannConstant = 3.1668114e-6; // E_h / K

MDIntegrator::MDIntegrator() : targetTemperature_(298.15 * kelvinToAtomicUnits) {
  setTimeStepInFemtoseconds(1.0);
  setRelaxationTimeFactor(10.0);
}

void MDIntegrator::setElementTypes(const Utils::ElementTypeCollection& elements) {
  nAtoms_ = static_cast<int>(elements.size());
  masses_ = Utils::Geometry::getMasses(elements);
  resetVelocities();
  resetAccelerations();
}

void MDIntegrator::calculateAccelerationsFromGradients(const Utils::GradientCollection& gradients) {
  assert(gradients.size() == accelerations_.size());
  for (int i = 0; i < nAtoms_; ++i)
    accelerations_.row(i) = (-1.0 / masses_[i]) * gradients.row(i);
}

void MDIntegrator::rescaleVelocitiesForTemperatureBath() {
  if (!temperatureBath_)
    return;
  double currentTemperature = getCurrentTemperature();
  double factor = std::sqrt(1. + timeStep_ / temperatureRelaxationTime_ * (targetTemperature_ / currentTemperature - 1));
  velocities_ *= factor;
}

void MDIntegrator::resetVelocities() {
  velocities_.resize(nAtoms_, 3);
  velocities_.setZero();
  sampleVelocitiesFromBoltzmannDistribution();
}

void MDIntegrator::resetAccelerations() {
  accelerations_.resize(nAtoms_, 3);
  accelerations_.setZero();
}

void MDIntegrator::setVelocities(const Utils::DisplacementCollection& velocities) {
  velocities_ = velocities;
}

Utils::DisplacementCollection MDIntegrator::getVelocities() const {
  return velocities_;
}

void MDIntegrator::setTimeStepInFemtoseconds(double fs) {
  timeStep_ = fs * femtosecondsToAtomicUnits / massConversionFactor;
  temperatureRelaxationTime_ = timeStep_ * relaxationTimeFactor_;
}

void MDIntegrator::setRelaxationTimeFactor(double factor) {
  relaxationTimeFactor_ = factor;
  temperatureRelaxationTime_ = timeStep_ * relaxationTimeFactor_;
}

void MDIntegrator::setTargetTemperatureInKelvin(double T) {
  targetTemperature_ = T * kelvinToAtomicUnits;
  temperatureBath_ = true;
}

double MDIntegrator::getCurrentTemperature() {
  int numberDegreesOfFreedom = 3 * nAtoms_; // TODO: Subtract translation and remove center of mass motion
  double currentTemperature =
      (velocities_.rowwise().squaredNorm().array() * Eigen::Map<const Eigen::ArrayXd>(masses_.data(), masses_.size())).sum() /
      numberDegreesOfFreedom;
  return currentTemperature;
}

void MDIntegrator::sampleVelocitiesFromBoltzmannDistribution() {
  std::mt19937 gen(seed_);
  double numerator = std::sqrt(targetTemperature_ * boltzmannConstant * atomicUnitsToKelvin);
  int index = 0;
  for (double mass : masses_) {
    std::normal_distribution<> d(0., numerator * std::sqrt(1. / mass));
    velocities_.row(index) = Position{d(gen), d(gen), d(gen)};
    ++index;
  }
}

void MDIntegrator::setSeed(int seed) {
  seed_ = seed;
}

void MDIntegrator::removeCenterOfMassLinearMomentum(const Eigen::MatrixX3d& positions) {
  const Eigen::VectorXd& masses = Eigen::Map<const Eigen::VectorXd>(masses_.data(), masses_.size());
  double totalMass = masses.sum();
  auto centerOfMass = Utils::Geometry::getCenterOfMass(positions, masses_);

  // Remove total linear momentum
  Eigen::MatrixX3d linearMomentumVector = velocities_.array().colwise() * masses.array();
  Eigen::RowVector3d centerOfMassVelocity = linearMomentumVector.colwise().sum() / totalMass;

  velocities_.rowwise() -= centerOfMassVelocity;
}

void MDIntegrator::removeCenterOfMassAngularMomentum(const Eigen::MatrixX3d& positions) {
  const Eigen::VectorXd& masses = Eigen::Map<const Eigen::VectorXd>(masses_.data(), masses_.size());
  double totalMass = masses.sum();
  auto centerOfMass = Utils::Geometry::getCenterOfMass(positions, masses_);
  auto inertiaTensor = Utils::Geometry::calculateInertiaTensor(positions, masses_, centerOfMass);
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