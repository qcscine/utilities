/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_DERIVATIVECOLLECTION_H
#define UTILSOS_DERIVATIVECOLLECTION_H

#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Utils/Math/AutomaticDifferentiation/Second3D.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @class DerivativeCollection DerivativeCollection.h
 * @brief This class is designed to store the first and optionally the second derivative of the
 *        energy with respect to the nuclear gradients.
 */
class DerivativeCollection {
 public:
  inline DerivativeCollection(int nAtoms, unsigned int order);
  /**
   * @brief Add another DerivativeCollection.
   */
  inline const DerivativeCollection& operator+=(const DerivativeCollection& dc);
  /**
   * @brief Reset all the derivatives to zero.
   */
  inline void setZero();
  /**
   * @brief Add a derivative object for the pair a->b.
   */
  inline void addDerivative(int a, int b, const AutomaticDifferentiation::Second3D& v);
  /**
   * @brief Getter for the Hessian Matrix.
   * @return HessianMatrix This ist just a type alias for Eigen::MatrixXd.
   */
  inline const HessianMatrix& getHessianMatrix() const;
  /**
   * @brief Return the reference gradients.
   * @return GradientCollection This ist just a type alias for Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>.
   */
  inline const GradientCollection& getReferenceGradients() const;
  /**
   * @brief Getter for th order of the derivative collection.
   * @return The derivative order.
   */
  inline unsigned int getOrder() const;
  /**
   * @brief Getter for the number of atoms.
   * @return The number of atoms.
   */
  inline int getNAtoms() const;

 private:
  // Hessian Matrix
  HessianMatrix hessianMatrix_;
  // Gradient Collection
  GradientCollection referenceGradients_;
  // Number of atoms
  int nAtoms_;
  // Derivative order.
  unsigned int order_;
};

inline DerivativeCollection::DerivativeCollection(int nAtoms, unsigned int order) : nAtoms_(nAtoms), order_(order) {
  if (order > 2) {
    throw std::runtime_error("Only derivatives up to order 2 are supported.");
  }
  referenceGradients_ = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>::Zero(nAtoms, 3);
  if (order > 1) {
    hessianMatrix_.resize(3 * nAtoms, 3 * nAtoms);
    hessianMatrix_.setZero();
  }
}

inline const DerivativeCollection& DerivativeCollection::operator+=(const DerivativeCollection& dc) {
  referenceGradients_ += dc.getReferenceGradients();
  if (order_ > 1 && dc.getOrder() > 1) {
    hessianMatrix_ += dc.getHessianMatrix();
  }
  return *this;
}

inline unsigned int DerivativeCollection::getOrder() const {
  return order_;
}

inline void DerivativeCollection::setZero() {
  referenceGradients_.setZero();
  hessianMatrix_.setZero();
}

inline void DerivativeCollection::addDerivative(int a, int b, const AutomaticDifferentiation::Second3D& v) {
  // Generate Eigen::Matrix3d from Second3D object:
  if (this->order_ > 1) {
    Eigen::Matrix3d m;
    m << v.XX(), v.XY(), v.XZ(), v.YX(), v.YY(), v.YZ(), v.ZX(), v.ZY(), v.ZZ();
    // Contributions added to the hessian:
    hessianMatrix_.block<3, 3>(3 * a, 3 * a) += m;
    hessianMatrix_.block<3, 3>(3 * b, 3 * b) += m;
    hessianMatrix_.block<3, 3>(3 * a, 3 * b) += -m;
    hessianMatrix_.block<3, 3>(3 * b, 3 * a) += -m;
  }
  // Contributions added to the gradients:
  referenceGradients_.row(b) += Gradient(v.dx(), v.dy(), v.dz());
  referenceGradients_.row(a) += Gradient(-v.dx(), -v.dy(), -v.dz());
}

inline const HessianMatrix& DerivativeCollection::getHessianMatrix() const {
  if (order_ < 2) {
    throw std::runtime_error("There is no Hessian available for derivatives collections with order 1 or lower.");
  }
  return hessianMatrix_;
}

inline const GradientCollection& DerivativeCollection::getReferenceGradients() const {
  return referenceGradients_;
}

inline int DerivativeCollection::getNAtoms() const {
  return nAtoms_;
}

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_DERIVATIVECOLLECTION_H
