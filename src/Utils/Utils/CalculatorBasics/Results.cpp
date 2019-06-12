/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/CalculatorBasics/Results.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include <boost/optional.hpp>

namespace Scine {
namespace Utils {

struct Results::Impl {
  boost::optional<std::string> description;
  boost::optional<double> energy;
  boost::optional<GradientCollection> gradients;
  boost::optional<HessianMatrix> hessian;
  boost::optional<Dipole> dipole;
  boost::optional<DipoleGradient> dipoleGradient;
  boost::optional<DipoleMatrix> dipoleMatrixAO, dipoleMatrixMO;
  boost::optional<Eigen::MatrixXd> oneElectronMatrix;
  boost::optional<SpinAdaptedMatrix> twoElectronMatrix;
  boost::optional<BondOrderCollection> bondOrders;
};

Results::Results() : pImpl_(std::make_unique<Impl>()) {
}

Results::Results(Results&& /*rhs*/) noexcept = default;

Results::Results(const Results& rhs) : Results() {
  if (rhs.hasDescription())
    setDescription(rhs.getDescription());
  if (rhs.hasEnergy())
    setEnergy(rhs.getEnergy());
  if (rhs.hasGradients())
    setGradients(rhs.getGradients());
  if (rhs.hasHessian())
    setHessian(rhs.getHessian());
  if (rhs.hasDipole())
    setDipole(rhs.getDipole());
  if (rhs.hasDipoleGradient())
    setDipoleGradient(rhs.getDipoleGradient());
  if (rhs.hasAODipoleMatrix())
    setAODipoleMatrix(rhs.getAODipoleMatrix());
  if (rhs.hasMODipoleMatrix())
    setMODipoleMatrix(rhs.getMODipoleMatrix());
  if (rhs.hasOneElectronMatrix())
    setOneElectronMatrix(rhs.getOneElectronMatrix());
  if (rhs.hasTwoElectronMatrix())
    setTwoElectronMatrix(rhs.getTwoElectronMatrix());
  if (rhs.hasBondOrders())
    setBondOrders(rhs.getBondOrders());
}

Results& Results::operator=(const Results& rhs) {
  static_assert(std::is_move_assignable<Results>::value, "Must be move assignable!");
  *this = Results(rhs);
  return *this;
}

Results::~Results() = default;

Results& Results::operator=(Results&& /*rhs*/) noexcept = default;

bool Results::hasDescription() const {
  return static_cast<bool>(pImpl_->description);
}

void Results::setDescription(std::string description) {
  pImpl_->description = std::move(description);
}

const std::string& Results::getDescription() const {
  if (!hasDescription()) {
    throw PropertyNotPresentException();
  }
  return *pImpl_->description;
}

bool Results::hasEnergy() const {
  return static_cast<bool>(pImpl_->energy);
}

void Results::setEnergy(double e) {
  pImpl_->energy = e;
}

double Results::getEnergy() const {
  if (!hasEnergy()) {
    throw PropertyNotPresentException();
  }
  return *pImpl_->energy;
}

bool Results::hasGradients() const {
  return static_cast<bool>(pImpl_->gradients);
}

void Results::setGradients(Utils::GradientCollection gradients) {
  pImpl_->gradients = std::move(gradients);
}

const Utils::GradientCollection& Results::getGradients() const {
  if (!hasGradients()) {
    throw PropertyNotPresentException();
  }
  return *pImpl_->gradients;
}

Utils::GradientCollection Results::takeGradients() {
  if (!hasGradients()) {
    throw PropertyNotPresentException();
  }
  auto gc = std::move(*pImpl_->gradients);
  pImpl_->gradients = boost::none;
  return gc;
}

bool Results::hasHessian() const {
  return static_cast<bool>(pImpl_->hessian);
}

void Results::setHessian(Utils::HessianMatrix hessian) {
  pImpl_->hessian = std::move(hessian);
}

const Utils::HessianMatrix& Results::getHessian() const {
  if (!hasHessian()) {
    throw PropertyNotPresentException();
  }
  return *pImpl_->hessian;
}

Utils::HessianMatrix Results::takeHessian() {
  if (!hasHessian()) {
    throw PropertyNotPresentException();
  }
  auto hessian = std::move(*pImpl_->hessian);
  pImpl_->hessian = boost::none;
  return hessian;
}

bool Results::hasDipole() const {
  return static_cast<bool>(pImpl_->dipole);
}

void Results::setDipole(Dipole dipole) {
  pImpl_->dipole = std::move(dipole);
};

const Dipole& Results::getDipole() const {
  if (!hasDipole()) {
    throw PropertyNotPresentException();
  }
  return *pImpl_->dipole;
}

bool Results::hasDipoleGradient() const {
  return static_cast<bool>(pImpl_->dipoleGradient);
}

void Results::setDipoleGradient(DipoleGradient dipoleGradient) {
  pImpl_->dipoleGradient = std::move(dipoleGradient);
}

const DipoleGradient& Results::getDipoleGradient() const {
  if (!hasDipoleGradient()) {
    throw PropertyNotPresentException();
  }
  return *pImpl_->dipoleGradient;
}

DipoleGradient Results::takeDipoleGradient() {
  if (!hasDipoleGradient()) {
    throw PropertyNotPresentException();
  }
  auto dipoleGradient = std::move(*pImpl_->dipoleGradient);
  pImpl_->dipoleGradient = boost::none;
  return dipoleGradient;
}

bool Results::hasAODipoleMatrix() const {
  return static_cast<bool>(pImpl_->dipoleMatrixAO);
}

void Results::setAODipoleMatrix(DipoleMatrix dipoleMatrixAO) {
  pImpl_->dipoleMatrixAO = std::move(dipoleMatrixAO);
}

const Utils::DipoleMatrix& Results::getAODipoleMatrix() const {
  if (!hasAODipoleMatrix()) {
    throw PropertyNotPresentException();
  }
  return *pImpl_->dipoleMatrixAO;
}

Utils::DipoleMatrix Results::takeAODipoleMatrix() {
  if (!hasAODipoleMatrix()) {
    throw PropertyNotPresentException();
  }
  auto dipoleMatrixAO = std::move(*pImpl_->dipoleMatrixAO);
  pImpl_->dipoleMatrixAO = boost::none;
  return dipoleMatrixAO;
}

bool Results::hasMODipoleMatrix() const {
  return static_cast<bool>(pImpl_->dipoleMatrixMO);
}

void Results::setMODipoleMatrix(DipoleMatrix dipoleMatrixMO) {
  pImpl_->dipoleMatrixMO = std::move(dipoleMatrixMO);
}

const Utils::DipoleMatrix& Results::getMODipoleMatrix() const {
  if (!hasMODipoleMatrix()) {
    throw PropertyNotPresentException();
  }
  return *pImpl_->dipoleMatrixMO;
}

Utils::DipoleMatrix Results::takeMODipoleMatrix() {
  if (!hasMODipoleMatrix()) {
    throw PropertyNotPresentException();
  }
  auto dipoleMatrixMO = std::move(*pImpl_->dipoleMatrixMO);
  pImpl_->dipoleMatrixMO = boost::none;
  return dipoleMatrixMO;
}

bool Results::hasOneElectronMatrix() const {
  return static_cast<bool>(pImpl_->oneElectronMatrix);
}

void Results::setOneElectronMatrix(Eigen::MatrixXd oneElectronMatrix) {
  pImpl_->oneElectronMatrix = std::move(oneElectronMatrix);
}

const Eigen::MatrixXd& Results::getOneElectronMatrix() const {
  if (!hasOneElectronMatrix()) {
    PropertyNotPresentException();
  }
  return *pImpl_->oneElectronMatrix;
}

Eigen::MatrixXd Results::takeOneElectronMatrix() {
  if (!hasOneElectronMatrix()) {
    throw PropertyNotPresentException();
  }
  auto oneElectronMatrix = std::move(*pImpl_->oneElectronMatrix);
  pImpl_->oneElectronMatrix = boost::none;
  return oneElectronMatrix;
}

bool Results::hasTwoElectronMatrix() const {
  return static_cast<bool>(pImpl_->twoElectronMatrix);
}

void Results::setTwoElectronMatrix(Utils::SpinAdaptedMatrix twoElectronMatrix) {
  pImpl_->twoElectronMatrix = std::move(twoElectronMatrix);
}

const Utils::SpinAdaptedMatrix& Results::getTwoElectronMatrix() const {
  if (!hasTwoElectronMatrix()) {
    throw PropertyNotPresentException();
  }
  return *pImpl_->twoElectronMatrix;
}

Utils::SpinAdaptedMatrix Results::takeTwoElectronMatrix() {
  if (!hasTwoElectronMatrix()) {
    throw PropertyNotPresentException();
  }
  auto twoElectronMatrix = std::move(*pImpl_->twoElectronMatrix);
  pImpl_->twoElectronMatrix = boost::none;
  return twoElectronMatrix;
}

bool Results::hasBondOrders() const {
  return static_cast<bool>(pImpl_->bondOrders);
}

void Results::setBondOrders(BondOrderCollection bondOrders) {
  pImpl_->bondOrders = std::move(bondOrders);
}

const BondOrderCollection& Results::getBondOrders() const {
  if (!hasBondOrders()) {
    throw PropertyNotPresentException();
  }
  return *pImpl_->bondOrders;
}

BondOrderCollection Results::takeBondOrders() {
  if (!hasBondOrders()) {
    throw PropertyNotPresentException();
  }
  auto bondOrders = std::move(*pImpl_->bondOrders);
  pImpl_->bondOrders = boost::none;
  return bondOrders;
}

} // namespace Utils
} // namespace Scine
