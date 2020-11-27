/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Geometry/InternalCoordinates.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Geometry/GeometryUtilities.h"
#include <libirc/irc.h>

namespace Scine {
namespace Utils {

struct InternalCoordinates::Impl {
  // External class
  std::unique_ptr<irc::IRC<Eigen::Vector3d, Eigen::VectorXd, Eigen::MatrixXd>> irc;
  // Backup transformation matrix
  std::unique_ptr<Eigen::MatrixXd> T;
};

InternalCoordinates::InternalCoordinates(const AtomCollection& atoms, bool rotTransOnly)
  : _oldCartesian(atoms.size() * 3) {
  irc::molecule::Molecule<Eigen::Vector3d> molecule;
  for (unsigned int i = 0; i < atoms.size(); i++) {
    auto symbol = ElementInfo::symbol(atoms.getElement(i));
    auto position = atoms.getPosition(i);
    molecule.push_back({symbol, {position[0], position[1], position[2]}});
    _oldCartesian[i * 3 + 0] = position[0];
    _oldCartesian[i * 3 + 1] = position[1];
    _oldCartesian[i * 3 + 2] = position[2];
  }
  this->_pImpl = std::make_unique<InternalCoordinates::Impl>();

  if (atoms.size() > 3 && !rotTransOnly) {
    try {
      this->_pImpl->irc = std::make_unique<irc::IRC<Eigen::Vector3d, Eigen::VectorXd, Eigen::MatrixXd>>(molecule);
      _oldInternal = this->_pImpl->irc->cartesian_to_irc(_oldCartesian);
    }
    catch (...) {
      // If not successful at least try to remove rotation and translation
      rotTransOnly = true;
    }
  }
  else {
    rotTransOnly = true;
  }
  if (rotTransOnly) {
    this->_pImpl->T = std::make_unique<Eigen::MatrixXd>(
        Geometry::calculateRotTransFreeTransformMatrix(atoms.getPositions(), atoms.getElements()));
  }
}

InternalCoordinates::~InternalCoordinates() = default;

PositionCollection InternalCoordinates::coordinatesToCartesian(const Eigen::VectorXd& internals, unsigned int maxIters,
                                                               double tolerance) const {
  if (this->_pImpl->T) {
    auto backtransformed = ((*this->_pImpl->T) * internals).eval();
    assert(static_cast<int>(backtransformed.size()) % 3 == 0);
    return Eigen::Map<PositionCollection>(backtransformed.data(), static_cast<int>(backtransformed.size() / 3), 3);
  }
  Eigen::VectorXd deltaInternal = (internals - _oldInternal).eval();
  auto result = this->_pImpl->irc->irc_to_cartesian(_oldInternal, deltaInternal, _oldCartesian, maxIters, tolerance);
  if (!result.converged)
    throw InternalCoordinatesException();
  _oldCartesian = result.x_c;
  _oldInternal = internals.eval();
  return Eigen::Map<PositionCollection>(_oldCartesian.data(), static_cast<int>(_oldCartesian.size() / 3), 3);
}

Eigen::VectorXd InternalCoordinates::coordinatesToInternal(const PositionCollection& cartesian) const {
  auto vector = Eigen::Map<const Eigen::VectorXd>(cartesian.data(), cartesian.size());
  if (this->_pImpl->T) {
    return this->_pImpl->T->transpose() * vector;
  }
  return this->_pImpl->irc->cartesian_to_irc(vector);
}

Eigen::VectorXd InternalCoordinates::gradientsToInternal(const GradientCollection& cartesian) const {
  auto vector = Eigen::Map<const Eigen::VectorXd>(cartesian.data(), cartesian.size());
  if (this->_pImpl->T) {
    return this->_pImpl->T->transpose() * vector;
  }
  return this->_pImpl->irc->grad_cartesian_to_projected_irc(vector);
}

Eigen::MatrixXd InternalCoordinates::inverseHessianGuess() const {
  if (this->_pImpl->T) {
    const size_t n = this->_pImpl->T->cols();
    return Eigen::MatrixXd::Identity(n, n);
  }
  return this->_pImpl->irc->projected_initial_hessian_inv();
}

Eigen::MatrixXd InternalCoordinates::hessianGuess() const {
  if (this->_pImpl->T) {
    const size_t n = this->_pImpl->T->cols();
    return Eigen::MatrixXd::Identity(n, n);
  }
  return this->_pImpl->irc->projected_initial_hessian();
}

Eigen::MatrixXd InternalCoordinates::projectHessianInverse(const Eigen::MatrixXd& internal) const {
  if (this->_pImpl->T) {
    return internal;
  }
  return this->_pImpl->irc->projected_hessian_inv(internal);
}

Eigen::MatrixXd InternalCoordinates::projectHessian(const Eigen::MatrixXd& internal) const {
  if (this->_pImpl->T) {
    return internal;
  }
  return this->_pImpl->irc->projected_hessian(internal);
}

} /* namespace Utils */
} /* namespace Scine */
