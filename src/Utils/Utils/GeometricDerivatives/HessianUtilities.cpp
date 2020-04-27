/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/GeometricDerivatives/HessianUtilities.h>
#include <Utils/Geometry/GeometryUtilities.h>
#include <Eigen/Dense>

namespace Scine {
namespace Utils {

HessianUtilities::HessianUtilities(const Eigen::MatrixXd& hessian, const ElementTypeCollection& elements,
                                   const PositionCollection& positions, bool massWeighted)
  : _hessian(hessian), _elements(elements), _positions(positions), _massWeighted(massWeighted) {
  _transformation = Geometry::calculateRotTransFreeTransformMatrix(positions, elements, _massWeighted);
}

void HessianUtilities::hessianUpdate() {
  _internalEValues.reset(nullptr);
  _internalEVectors.reset(nullptr);
}

void HessianUtilities::hessianUpdate(const HessianMatrix& hessian) {
  _hessian = std::ref(hessian);
  this->hessianUpdate();
}

const Eigen::MatrixXd& HessianUtilities::getTransformationMatrix() const {
  return _transformation;
}

const Eigen::VectorXd& HessianUtilities::getInternalEigenvalues() {
  if (!_internalEValues)
    this->calculateInternal();
  return *_internalEValues;
}

const Eigen::MatrixXd& HessianUtilities::getInternalEigenvectors() {
  if (!_internalEVectors)
    this->calculateInternal();
  return *_internalEVectors;
}

Eigen::MatrixXd HessianUtilities::getBackTransformedInternalEigenvectors() {
  if (!_internalEVectors) {
    this->calculateInternal();
  }
  if (_massWeighted) {
    auto masses = Geometry::getMasses(_elements);
    Eigen::MatrixXd ret = _transformation * (*_internalEVectors);
    for (int i = 0; i < masses.size(); ++i)
      ret.middleRows(3 * i, 3) *= 1. / std::sqrt(masses[i]);
    ret.colwise().normalize();
    return ret;
  }
  else {
    return _transformation * (*_internalEVectors);
  }
}

Eigen::MatrixXd HessianUtilities::getInternalHessian() const {
  if (_massWeighted) {
    Eigen::MatrixXd m = _hessian.get();
    auto masses = Geometry::getMasses(_elements);
    for (unsigned i = 0; i < masses.size(); ++i) {
      double invSqrtM = 1. / std::sqrt(masses[i]);
      m.middleCols(3 * i, 3) *= invSqrtM;
      m.middleRows(3 * i, 3) *= invSqrtM;
    }
    return _transformation.transpose() * m * _transformation;
  }
  else {
    return _transformation.transpose() * _hessian.get() * _transformation;
  }
}

void HessianUtilities::calculateInternal() {
  if (_massWeighted) {
    Eigen::MatrixXd m = _hessian.get();
    auto masses = Geometry::getMasses(_elements);
    // Get mass weighted Hessian
    for (unsigned i = 0; i < masses.size(); ++i) {
      double invSqrtM = 1. / std::sqrt(masses[i]);
      m.middleCols(3 * i, 3) *= invSqrtM;
      m.middleRows(3 * i, 3) *= invSqrtM;
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    // Project out by multiplying by  transformation matrix (has to be mass-weighted, too!)
    es.compute(_transformation.transpose() * m * _transformation);
    _internalEValues = std::make_unique<Eigen::VectorXd>(es.eigenvalues());
    _internalEVectors = std::make_unique<Eigen::MatrixXd>(es.eigenvectors());
  }
  else {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(_transformation.transpose() * _hessian.get() * _transformation);
    _internalEValues = std::make_unique<Eigen::VectorXd>(es.eigenvalues());
    _internalEVectors = std::make_unique<Eigen::MatrixXd>(es.eigenvectors());
  }
}

} // namespace Utils
} // namespace Scine
