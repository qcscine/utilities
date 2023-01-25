/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_SPINADAPTEDVECTOR_H
#define UTILSOS_SPINADAPTEDVECTOR_H

#include <Eigen/Core>
#include <utility>

namespace Scine {
namespace Utils {

/*!
 * Class defining a matrix and which can be used for both spin-restricted and spin-unrestricted formalisms in electronic
 * structure calculation methods. There is only slight overhead if only the restricted formulation is needed.
 */

class SpinAdaptedVector {
 public:
  using Vector = Eigen::VectorXd;

  /*! Construct a SpinAdaptedVector with restricted part only from an Eigen::VectorXd. */
  static SpinAdaptedVector createRestricted(Eigen::VectorXd m);
  /*! Construct a SpinAdaptedVector with alpha and beta parts only from two Eigen::VectorXd. */
  static SpinAdaptedVector createUnrestricted(Eigen::VectorXd alpha, Eigen::VectorXd beta);

  template<typename T>
  void setRestrictedVector(T&& restrictedVector);
  template<typename T>
  void setAlphaVector(T&& alphaVector);
  template<typename T>
  void setBetaVector(T&& betaVector);

  /*! \name Accessors to the underlying matrices by reference.  @{ */
  Vector& restrictedVector();
  const Vector& restrictedVector() const;
  Vector& alphaVector();
  const Vector& alphaVector() const;
  Vector& betaVector();
  const Vector& betaVector() const;
  /*! @} */

  /*! \name Accessors to the matrix elements.
      Using template functions to allow perfect forwarding to the Eigen functions.  @{ */
  template<typename Index>
  double restricted(Index i) const;
  template<typename Index>
  double alpha(Index i) const;
  template<typename Index>
  double beta(Index i) const;
  /*! @} */
  /*! \name Setters for the matrix elements.
      Using template functions to allow perfect forwarding to the Eigen functions.  @{ */
  template<typename Index>
  void setRestricted(Index i, double x);
  template<typename Index>
  void setAlpha(Index i, double x);
  template<typename Index>
  void setBeta(Index i, double x);
  /*! @} */

  /*! Set the size of the matrices (number of atomic orbitals) */
  void resize(int nAOs);

 protected:
  Vector restrictedVector_;
  Vector alphaVector_;
  Vector betaVector_;
};

template<typename T>
void SpinAdaptedVector::setRestrictedVector(T&& restrictedVector) {
  restrictedVector_ = std::forward<T>(restrictedVector);
}

template<typename T>
void SpinAdaptedVector::setAlphaVector(T&& alphaVector) {
  alphaVector_ = std::forward<T>(alphaVector);
}

template<typename T>
void SpinAdaptedVector::setBetaVector(T&& betaVector) {
  betaVector_ = std::forward<T>(betaVector);
}

template<typename Index>
double SpinAdaptedVector::restricted(Index i) const {
  return restrictedVector_(i);
}

template<typename Index>
double SpinAdaptedVector::alpha(Index i) const {
  return alphaVector_(i);
}

template<typename Index>
double SpinAdaptedVector::beta(Index i) const {
  return betaVector_(i);
}

template<typename Index>
void SpinAdaptedVector::setRestricted(Index i, double x) {
  restrictedVector_(i) = x;
}

template<typename Index>
void SpinAdaptedVector::setAlpha(Index i, double x) {
  alphaVector_(i) = x;
}

template<typename Index>
void SpinAdaptedVector::setBeta(Index i, double x) {
  betaVector_(i) = x;
}

inline SpinAdaptedVector::Vector& SpinAdaptedVector::restrictedVector() {
  return restrictedVector_;
}

inline const SpinAdaptedVector::Vector& SpinAdaptedVector::restrictedVector() const {
  return restrictedVector_;
}

inline SpinAdaptedVector::Vector& SpinAdaptedVector::alphaVector() {
  return alphaVector_;
}

inline const SpinAdaptedVector::Vector& SpinAdaptedVector::alphaVector() const {
  return alphaVector_;
}

inline SpinAdaptedVector::Vector& SpinAdaptedVector::betaVector() {
  return betaVector_;
}

inline const SpinAdaptedVector::Vector& SpinAdaptedVector::betaVector() const {
  return betaVector_;
}

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_SPINADAPTEDVECTOR_H
