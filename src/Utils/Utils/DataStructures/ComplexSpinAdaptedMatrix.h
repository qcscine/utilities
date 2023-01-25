/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_COMPLEXSPINADAPTEDMATRIX_H
#define UTILS_COMPLEXSPINADAPTEDMATRIX_H

#include <Eigen/Core>
#include <utility>

namespace Scine {
namespace Utils {

/*!
 * Class defining a matrix and which can be used for both spin-restricted and spin-unrestricted formalisms in electronic
 * structure calculation methods. There is only slight overhead if only the restricted formulation is needed.
 */

class ComplexSpinAdaptedMatrix {
 public:
  using Matrix = Eigen::MatrixXcd;

  /*! Construct a ComplexSpinAdaptedMatrix with restricted part only from an Eigen::MatrixXd. */
  static ComplexSpinAdaptedMatrix createRestricted(Eigen::MatrixXd m);
  /*! Construct a ComplexSpinAdaptedMatrix with alpha and beta parts only from two Eigen::MatrixXd. */
  static ComplexSpinAdaptedMatrix createUnrestricted(Eigen::MatrixXd alpha, Eigen::MatrixXd beta);

  template<typename T>
  void setRestrictedMatrix(T&& restrictedMatrix);
  template<typename T>
  void setAlphaMatrix(T&& alphaMatrix);
  template<typename T>
  void setBetaMatrix(T&& betaMatrix);

  /*! \name Accessors to the underlying matrices by reference.  @{ */
  Matrix& restrictedMatrix();
  const Matrix& restrictedMatrix() const;
  Matrix& alphaMatrix();
  const Matrix& alphaMatrix() const;
  Matrix& betaMatrix();
  const Matrix& betaMatrix() const;
  /*! @} */

  /*! \name Accessors to the matrix elements.
      Using template functions to allow perfect forwarding to the Eigen functions.  @{ */
  template<typename Index>
  double restricted(Index i, Index j) const;
  template<typename Index>
  double alpha(Index i, Index j) const;
  template<typename Index>
  double beta(Index i, Index j) const;
  /*! @} */
  /*! \name Setters for the matrix elements.
      Using template functions to allow perfect forwarding to the Eigen functions.  @{ */
  template<typename Index>
  void setRestricted(Index i, Index j, double x);
  template<typename Index>
  void setAlpha(Index i, Index j, double x);
  template<typename Index>
  void setBeta(Index i, Index j, double x);
  /*! @} */

  /*! Set the size of the matrices (number of atomic orbitals) */
  void resize(int nAOs);

 protected:
  Matrix restrictedMatrix_;
  Matrix alphaMatrix_;
  Matrix betaMatrix_;
};

template<typename T>
void ComplexSpinAdaptedMatrix::setRestrictedMatrix(T&& restrictedMatrix) {
  restrictedMatrix_ = std::forward<T>(restrictedMatrix);
}

template<typename T>
void ComplexSpinAdaptedMatrix::setAlphaMatrix(T&& alphaMatrix) {
  alphaMatrix_ = std::forward<T>(alphaMatrix);
}

template<typename T>
void ComplexSpinAdaptedMatrix::setBetaMatrix(T&& betaMatrix) {
  betaMatrix_ = std::forward<T>(betaMatrix);
}

template<typename Index>
double ComplexSpinAdaptedMatrix::restricted(Index i, Index j) const {
  return restrictedMatrix_(i, j);
}

template<typename Index>
double ComplexSpinAdaptedMatrix::alpha(Index i, Index j) const {
  return alphaMatrix_(i, j);
}

template<typename Index>
double ComplexSpinAdaptedMatrix::beta(Index i, Index j) const {
  return betaMatrix_(i, j);
}

template<typename Index>
void ComplexSpinAdaptedMatrix::setRestricted(Index i, Index j, double x) {
  restrictedMatrix_(i, j) = x;
}

template<typename Index>
void ComplexSpinAdaptedMatrix::setAlpha(Index i, Index j, double x) {
  alphaMatrix_(i, j) = x;
}

template<typename Index>
void ComplexSpinAdaptedMatrix::setBeta(Index i, Index j, double x) {
  betaMatrix_(i, j) = x;
}

inline ComplexSpinAdaptedMatrix::Matrix& ComplexSpinAdaptedMatrix::restrictedMatrix() {
  return restrictedMatrix_;
}

inline const ComplexSpinAdaptedMatrix::Matrix& ComplexSpinAdaptedMatrix::restrictedMatrix() const {
  return restrictedMatrix_;
}

inline ComplexSpinAdaptedMatrix::Matrix& ComplexSpinAdaptedMatrix::alphaMatrix() {
  return alphaMatrix_;
}

inline const ComplexSpinAdaptedMatrix::Matrix& ComplexSpinAdaptedMatrix::alphaMatrix() const {
  return alphaMatrix_;
}

inline ComplexSpinAdaptedMatrix::Matrix& ComplexSpinAdaptedMatrix::betaMatrix() {
  return betaMatrix_;
}

inline const ComplexSpinAdaptedMatrix::Matrix& ComplexSpinAdaptedMatrix::betaMatrix() const {
  return betaMatrix_;
}

} // namespace Utils
} // namespace Scine
#endif // UTILS_SPINADAPTEDMATRIX_H
