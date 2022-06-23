/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/GeometricDerivatives/HessianUtilities.h>
#include <Utils/Geometry/GeometryUtilities.h>
#include <Eigen/Dense>

namespace Scine {
namespace Utils {

HessianUtilities::HessianUtilities(const Eigen::MatrixXd& hessian, const ElementTypeCollection& elements,
                                   const PositionCollection& positions, bool massWeighted)
  : _massWeighted(massWeighted), _hessian(hessian), _elements(elements) {
  _transformation = Geometry::Transformations::calculateRotTransFreeTransformMatrix(positions, elements, _massWeighted);
}

HessianUtilities::HessianUtilities(const Eigen::MatrixXd& hessian, const ElementTypeCollection& elements,
                                   const PositionCollection& positions, const GradientCollection& gradient, bool massWeighted)
  : _massWeighted(massWeighted), _hessian(hessian), _elements(elements) {
  _transformation =
      Geometry::Transformations::calculateRotTransFreeTransformMatrix(positions, elements, gradient, _massWeighted);
}

void HessianUtilities::hessianUpdate() {
  _internalEValues.reset(nullptr);
  _internalEVectors.reset(nullptr);
  _gradient.reset(nullptr);
}

void HessianUtilities::hessianUpdate(const HessianMatrix& hessian) {
  _hessian = std::ref(hessian);
  this->hessianUpdate();
}

const Eigen::MatrixXd& HessianUtilities::getTransformationMatrix() const {
  return _transformation;
}

const Eigen::VectorXd& HessianUtilities::getInternalEigenvalues() {
  if (!_internalEValues) {
    this->calculateInternal();
  }
  return *_internalEValues;
}

const Eigen::MatrixXd& HessianUtilities::getInternalEigenvectors() {
  if (!_internalEVectors) {
    this->calculateInternal();
  }
  return *_internalEVectors;
}

Eigen::MatrixXd HessianUtilities::getBackTransformedInternalEigenvectors(bool normalize) {
  if (!_internalEVectors) {
    this->calculateInternal();
  }
  if (_massWeighted) {
    auto masses = Geometry::Properties::getMasses(_elements);
    Eigen::MatrixXd ret = _transformation * (*_internalEVectors);
    // Back-scale the mass-weighted coordinates to cartesian coordinates.
    const int N = masses.size();
    for (int i = 0; i < N; ++i) {
      ret.middleRows(3 * i, 3) *= 1. / std::sqrt(masses[i]);
    }
    if (normalize)
      ret.colwise().normalize();
    return ret;
  }

  return _transformation * (*_internalEVectors);
}

Eigen::MatrixXd HessianUtilities::getInternalHessian() const {
  if (_massWeighted) {
    Eigen::MatrixXd m = _hessian.get();
    auto masses = Geometry::Properties::getMasses(_elements);
    for (unsigned i = 0; i < masses.size(); ++i) {
      double invSqrtM = 1. / std::sqrt(masses[i]);
      m.middleCols(3 * i, 3) *= invSqrtM;
      m.middleRows(3 * i, 3) *= invSqrtM;
    }
    return _transformation.transpose() * m * _transformation;
  }

  return _transformation.transpose() * _hessian.get() * _transformation;
}

void HessianUtilities::calculateInternal() {
  // Check whether internal Hessian would be empty
  if (_transformation.size() == 0) {
    throw EmptyInternalHessianException();
  }

  if (_massWeighted) {
    Eigen::MatrixXd m = _hessian.get();
    auto masses = Geometry::Properties::getMasses(_elements);
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
