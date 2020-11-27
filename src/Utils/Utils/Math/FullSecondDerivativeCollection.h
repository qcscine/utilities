/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_FULLSECONDDERIVATIVECOLLECTION_H
#define UTILS_MATH_FULLSECONDDERIVATIVECOLLECTION_H

#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Utils/Math/AutomaticDifferentiation/Second3D.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @class FullSecondDerivativeCollection FullSecondDerivativeCollection.h
 * @brief Container class for full second derivatives (Hessian matrix + first derivatives).
 */
class FullSecondDerivativeCollection {
 public:
  /**
   * @brief Constructor.
   * @param N Size of the container (number of atoms).
   */
  explicit FullSecondDerivativeCollection(int N = 0);

  /**
   * @brief Add another FullSecondDerivativeCollection.
   */
  const FullSecondDerivativeCollection& operator+=(const FullSecondDerivativeCollection& dc);
  /**
   * @brief Reset all the derivatives to zero.
   */
  void setZero();
  /**
   * @brief Generate a GradientCollection from the full second derivatives and the displacements.
   */
  GradientCollection generateGradients(const DisplacementCollection& displacements) const;
  /**
   * @brief Generate a Gradient from the full second derivatives and the displacements.
   */
  Gradient generateGradient(int idx, const DisplacementCollection& displacements) const;
  /**
   * @brief Add a derivative object for the pair a->b.
   */
  void addDerivative(int a, int b, const AutomaticDifferentiation::Second3D& v);
  /**
   * @brief Getter for the Hessian Matrix.
   * @return HessianMatrix This ist just a type alias for Eigen::MatrixXd.
   */
  const HessianMatrix& getHessianMatrix() const;
  /**
   * @brief Return the reference gradients.
   * @return GradientCollection This ist just a type alias for Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>.
   */
  const GradientCollection& getReferenceGradients() const;

 private:
  // Getter for the hessian that returns a reference, not a const reference.
  Eigen::MatrixXd& hessian();

  // Hessian Matrix
  HessianMatrix hessianMatrix_;
  // Gradient Collection
  GradientCollection referenceGradients_;
  // Number of atoms
  int nAtoms_;
};

// Constructor with the number of atoms given as an argument
inline FullSecondDerivativeCollection::FullSecondDerivativeCollection(int N) {
  // Set number of atoms
  nAtoms_ = N;
  // Initiate hessian matrix
  hessian() = Eigen::MatrixXd::Zero(3 * N, 3 * N);
  // Initiate gradients
  referenceGradients_.resize(N, 3);
}

// Operator for adding two containers
inline const FullSecondDerivativeCollection& FullSecondDerivativeCollection::operator+=(const FullSecondDerivativeCollection& dc) {
  hessian() += dc.getHessianMatrix();
  referenceGradients_ += dc.referenceGradients_;
  return *this;
}

// Set all derivatives zero
inline void FullSecondDerivativeCollection::setZero() {
  hessian().setZero();
  referenceGradients_.setZero();
}

// This function generates gradients from displacements
inline GradientCollection FullSecondDerivativeCollection::generateGradients(const DisplacementCollection& displacements) const {
  GradientCollection gc(nAtoms_, 3);
  for (int i = 0; i < nAtoms_; ++i) {
    Eigen::Vector3d v = referenceGradients_.row(i);
    for (int j = 0; j < nAtoms_; ++j)
      v += getHessianMatrix().block<3, 3>(3 * i, 3 * j) * displacements.row(j).transpose();
    gc.row(i) = v;
  }

  return gc;
}

inline Gradient FullSecondDerivativeCollection::generateGradient(int /*idx*/, const DisplacementCollection& /*unused*/) const {
  throw std::runtime_error("Not implemented yet");
}

// This functions allows for adding a Second3D object that is refering to the atoms a and b to the container.
inline void FullSecondDerivativeCollection::addDerivative(int a, int b, const AutomaticDifferentiation::Second3D& v) {
  // Generate Eigen::Matrix3d from Second3D object:
  Eigen::Matrix3d m;
  m << v.XX(), v.XY(), v.XZ(), v.YX(), v.YY(), v.YZ(), v.ZX(), v.ZY(), v.ZZ();
  // Contributions added to the hessian:
  hessian().block<3, 3>(3 * a, 3 * a) += m;
  hessian().block<3, 3>(3 * b, 3 * b) += m;
  hessian().block<3, 3>(3 * a, 3 * b) += -m;
  hessian().block<3, 3>(3 * b, 3 * a) += -m;
  // Contributions added to the gradients:
  referenceGradients_.row(b) += Gradient(v.dx(), v.dy(), v.dz());
  referenceGradients_.row(a) += Gradient(-v.dx(), -v.dy(), -v.dz());
}

// Getter for hessian matrix (by const reference)
inline const HessianMatrix& FullSecondDerivativeCollection::getHessianMatrix() const {
  return hessianMatrix_;
}

// Getter for hessian matrix (by reference)
inline Eigen::MatrixXd& FullSecondDerivativeCollection::hessian() {
  return hessianMatrix_;
}

// Getter for gradients
inline const GradientCollection& FullSecondDerivativeCollection::getReferenceGradients() const {
  return referenceGradients_;
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_FULLSECONDDERIVATIVECOLLECTION_H
