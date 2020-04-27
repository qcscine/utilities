/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_DENSITYMATRIX_H
#define UTILS_DENSITYMATRIX_H

#include <Utils/DataStructures/SpinAdaptedMatrix.h>

namespace Scine {
namespace Utils {

/*!
 * Class defining a density matrix for use in electronic structure calculation methods.
 * Is adequate for both restricted and unrestricted formulations.
 * Allows for an easy transfer of density matrices between different instances of a quantum chemical method.
 * There is no overhead if only the restricted formulation is needed.
 */
class DensityMatrix {
 public:
  using Matrix = Eigen::MatrixXd;

  DensityMatrix() = default;

  /*! Set the restricted density matrix.
      If in unrestricted mode, the alpha and beta density are set to half each. */
  void setDensity(Matrix&& restrictedMatrix, int nElectrons);
  /*! Set the unrestricted density matrices.
      The "total" density matrix is automatically updated. */
  void setDensity(Matrix&& alphaMatrix, Matrix&& betaMatrix, int nAlphaElectrons, int nBetaElectrons);

  /*! \name Accessors to the underlying matrices by const reference.  @{ */
  const Matrix& restrictedMatrix() const;
  const Matrix& alphaMatrix() const;
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

  /*! Set whether the unrestricted formalism is to use. */
  void setUnrestricted(bool b);
  /*! Get whether the unrestricted formalism is used. */
  bool unrestricted() const;
  /*! Get whether the restricted formalism is used. */
  bool restricted() const;
  /*! Resizes the alpha and beta density matrices and sets them each to half the full density matrix. */
  void setAlphaAndBetaFromRestrictedDensity();

  /*! Set the size of the matrices (number of atomic orbitals) */
  void resize(int nAOs);

  /*! Returns the total number of electrons in the density matrix.
   * Integer version of getOccupation().*/
  int numberElectrons() const;
  /*! Return the number of electrons corresponding to the alpha density matrix.
   * Integer version of getAlphaOccupation(). */
  int numberElectronsInAlphaMatrix() const;
  /*! Return the number of electrons corresponding to the beta density matrix.
   * Integer version of getBetaOccupation().*/
  int numberElectronsInBetaMatrix() const;
  double getOccupation() const;
  double getAlphaOccupation() const;
  double getBetaOccupation() const;

  /*! \name Arithmetic operations
   * In the case of functions for pairs of density matrices, the precondition is that they both have the same size
   * and are either RHF-RHF or UHF-UHF.
   * @{ */
  DensityMatrix& operator+=(const DensityMatrix& rhs);
  DensityMatrix operator+(const DensityMatrix& rhs) const;
  DensityMatrix& operator-=(const DensityMatrix& rhs);
  DensityMatrix operator-(const DensityMatrix& rhs) const;
  DensityMatrix& operator*=(double f);
  DensityMatrix operator*(double f);
  /*! @} */

 private:
  SpinAdaptedMatrix matrix_;
  bool unrestricted_ = false; /*!< true if the alpha and beta densities are needed and were set. */
  double alphaOccupation_ = 0.0, betaOccupation_ = 0.0;
};

inline bool DensityMatrix::unrestricted() const {
  return unrestricted_;
}

inline bool DensityMatrix::restricted() const {
  return !unrestricted();
}

template<typename Index>
double DensityMatrix::restricted(Index i, Index j) const {
  return restrictedMatrix()(i, j);
}

template<typename Index>
double DensityMatrix::alpha(Index i, Index j) const {
  return alphaMatrix()(i, j);
}

template<typename Index>
double DensityMatrix::beta(Index i, Index j) const {
  return betaMatrix()(i, j);
}

inline const DensityMatrix::Matrix& DensityMatrix::restrictedMatrix() const {
  return matrix_.restrictedMatrix();
}

inline const DensityMatrix::Matrix& DensityMatrix::alphaMatrix() const {
  return matrix_.alphaMatrix();
}

inline const DensityMatrix::Matrix& DensityMatrix::betaMatrix() const {
  return matrix_.betaMatrix();
}

inline int DensityMatrix::numberElectrons() const {
  return static_cast<int>(std::lround(getOccupation()));
}

inline int DensityMatrix::numberElectronsInAlphaMatrix() const {
  return static_cast<int>(std::lround(getAlphaOccupation()));
}

inline int DensityMatrix::numberElectronsInBetaMatrix() const {
  return static_cast<int>(std::lround(getBetaOccupation()));
}

inline double DensityMatrix::getOccupation() const {
  return getAlphaOccupation() + getBetaOccupation();
}

inline double DensityMatrix::getAlphaOccupation() const {
  return alphaOccupation_;
}

inline double DensityMatrix::getBetaOccupation() const {
  return betaOccupation_;
}

} // namespace Utils
} // namespace Scine
#endif // UTILS_DENSITYMATRIX_H