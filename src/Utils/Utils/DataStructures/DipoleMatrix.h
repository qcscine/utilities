/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_DIPOLEMATRIX_H
#define UTILS_DIPOLEMATRIX_H

#include "MatrixWithDerivatives.h"
#include <Eigen/Core>
#include <array>

namespace Scine {
namespace Utils {

class InvalidDimensionSizeException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "Invalid dimension index. Possible dimensions: 0 = x, 1 = y, 2 = z for dipole matrix"
           " and 0 = xx, 1 = xy, 2 = xz, 3 = yy, 4 = yz, 5 = zz for quadrupole matrix.";
  }
};

/**
 * @class DipoleMatrix @file DipoleMatrix.h
 * @brief Class representing dipole matrix integrals with derivatives up to 2nd order.
 * This class stores the dipole integrals \f$\langle mu|r|nu \rangle\f$ and their derivatives in their three components
 x, y, z.
 *
 * The derivatives are accessible through the .derivatives() method of the
 * underlying AutomaticDifferentiation::First3D/Second3D class.
 *
 * @code{cpp}
   int arbitraryDimension = 10;
   DipoleMatrix dippoleMatrixZero, dipoleMatrixFirst;

   // Fill DipoleMatrix with random values and derivatives.
   using DerivType = Eigen::Matrix<AutomaticDifferentiation::First3D, Eigen::Dynamic, Eigen::Dynamic>;

   Eigen::MatrixXd XMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   Eigen::MatrixXd YMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   Eigen::MatrixXd ZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   DerivType Xderivatives, Yderivatives, Zderivatives;
   Xderivatives.resize(arbitraryDimension, arbitraryDimension);
   Yderivatives.resize(arbitraryDimension, arbitraryDimension);
   Zderivatives.resize(arbitraryDimension, arbitraryDimension);
   for (int i = 0; i < arbitraryDimension; ++i) {
     for (int j = 0; j < arbitraryDimension; ++j) {
       Xderivatives(i, j) = {XMatrix(i, j), Eigen::Vector3d::Random()};
       Yderivatives(i, j) = {YMatrix(i, j), Eigen::Vector3d::Random()};
       Zderivatives(i, j) = {ZMatrix(i, j), Eigen::Vector3d::Random()};
     }
   }
   // Fill zeroth order dipole matrix
   dipoleMatrixZero.x().get<DerivativeOrder::Zero>(XMatrix);
   dipoleMatrixZero.y().get<DerivativeOrder::Zero>(YMatrix);
   dipoleMatrixZero.z().get<DerivativeOrder::Zero>(ZMatrix);

   // Fill first order dipole matrix with the first derivative integrals
   dipoleMatrixFirst.x().get<DerivativeOrder::One>(Xderivatives);
   dipoleMatrixFirst.y().get<DerivativeOrder::One>(Yderivatives);
   dipoleMatrixFirst.z().get<DerivativeOrder::One>(Zderivatives);

   // Get the zeroth and first derivatives of the first element of the first derivative matrix
   dipoleMatrixZero.x().get<DerivativeOrder::Zero>()(0, 0);
   // Get the zeroth and first derivatives of the first element of the first derivative matrix
   AutomaticDifferentiation::getValue3DAsDouble(dipoleMatrixFirst.x().get<DerivativeOrder::One>()(0, 0));
   // or
   dipoleMatrixFirst.x().get<DerivativeOrder::One>()(0, 0).value();
   dipoleMatrixFirst.x().get<DerivativeOrder::One>()(0, 0).derivatives();
 * @endcode
 *
 */
class DipoleMatrix {
 public:
  /**
   * @brief Rule of 6
   * @{
   */
  DipoleMatrix();
  DipoleMatrix(const DipoleMatrix& rhs);
  DipoleMatrix& operator=(const DipoleMatrix& rhs);
  DipoleMatrix(DipoleMatrix&& rhs) noexcept;
  DipoleMatrix& operator=(DipoleMatrix&& rhs) noexcept;
  ~DipoleMatrix();
  //! @}
  /**
   * @brief Resets the dipole matrix.
   * @param dimension the number of basis function in the system, i.e. the dimension of the integral matrices.
   */
  void reset(int dimension);

  /**
   * @brief Setters for the x, y and z components of the matrix.
   * @{
   */
  MatrixWithDerivatives& x() {
    return dipoleMatrices_[0];
  }
  MatrixWithDerivatives& y() {
    return dipoleMatrices_[1];
  }
  MatrixWithDerivatives& z() {
    return dipoleMatrices_[2];
  }
  //! @}
  /**
   * @brief Getters for the x, y and z components of the matrix.
   * @{
   */
  const MatrixWithDerivatives& x() const {
    return dipoleMatrices_[0];
  }
  const MatrixWithDerivatives& y() const {
    return dipoleMatrices_[1];
  }
  const MatrixWithDerivatives& z() const {
    return dipoleMatrices_[2];
  }
  //! @}
  /**
   * @brief Returns the zeroth derivative components by index. 0 is x, 1 is y and 2 is z.
   * @{
   */
  Eigen::MatrixXd& operator[](int index);
  const Eigen::MatrixXd& operator[](int index) const;
  /**
   * @brief Returns the derivative components by index. 0 is x, 1 is y and 2 is z.
   * @{
   */
  MatrixWithDerivatives& at(int index) {
    return dipoleMatrices_.at(index);
  }

  const MatrixWithDerivatives& at(int index) const {
    return dipoleMatrices_.at(index);
  }
  //! @}
  // TODO: Implement symmetrization.
 private:
  std::array<MatrixWithDerivatives, 3> dipoleMatrices_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_DIPOLEMATRIX_H
