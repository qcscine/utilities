/**
 * @file QuadrupoleMatrix.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_QUADRUPOLEMATRIX_H
#define UTILSOS_QUADRUPOLEMATRIX_H

#include "MatrixWithDerivatives.h"
#include <array>

namespace Scine {
namespace Utils {

/**
 * @class QuadrupoleMatrix @file QuadrupoleMatrix.h
 * @brief Class representing quadrupole matrix integrals with derivatives up to 2nd order.
 *
 * This class stores the quadrupole integrals \f$\langle mu|q_{rr}|nu \rangle\f$
 * and their derivatives in their three components x, y, z.
 *
 * The class stores only the upper diagonal of the quadrupole moment tensor:
 * xx, xy, xz, yy, yz, zz.
 * The derivatives are accessible through the .derivatives() method of the
 * underlying AutomaticDifferentiation::First3D/Second3D class.
 * The templated getters and setters allow for efficient witing of
 * generic code, independent from the desired derivative order.
 * One just has to have specializations for all orders.
 *
 * @code{cpp}
   int arbitraryDimension = 10;
   QuadrupoleMatrix quadrupoleMatrixZero, quadrupoleMatrixFirst;

   // Fill QuadrupoleMatrix with random values and derivatives.
   DerivType = Eigen::Matrix<AutomaticDifferentiation::First3D, Eigen::Dynamic, Eigen::Dynamic>;

   Eigen::MatrixXd XXMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   Eigen::MatrixXd XYMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   Eigen::MatrixXd XZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   Eigen::MatrixXd YYMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   Eigen::MatrixXd YZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   Eigen::MatrixXd ZZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   DerivType XXderivatives, XYderivatives, XZderivatives,
             YYderivatives, YZderivatives, ZZderivatives;
   XXderivatives.resize(arbitraryDimension, arbitraryDimension);
   XYderivatives.resize(arbitraryDimension, arbitraryDimension);
   XZderivatives.resize(arbitraryDimension, arbitraryDimension);
   YYderivatives.resize(arbitraryDimension, arbitraryDimension);
   YZderivatives.resize(arbitraryDimension, arbitraryDimension);
   ZZderivatives.resize(arbitraryDimension, arbitraryDimension);
   for (int i = 0; i < arbitraryDimension; ++i) {
     for (int j = 0; j < arbitraryDimension; ++j) {
       XXderivatives(i, j) = {XXMatrix(i, j), Eigen::Vector3d::Random()};
       XYderivatives(i, j) = {XYMatrix(i, j), Eigen::Vector3d::Random()};
       XZderivatives(i, j) = {XZMatrix(i, j), Eigen::Vector3d::Random()};
       YYderivatives(i, j) = {YYMatrix(i, j), Eigen::Vector3d::Random()};
       YZderivatives(i, j) = {YZMatrix(i, j), Eigen::Vector3d::Random()};
       ZZderivatives(i, j) = {ZZMatrix(i, j), Eigen::Vector3d::Random()};
     }
   }
   // Fill with the zero-th derivative integrals
   quadrupoleMatrixZero.xx().get<DerivativeOrder::Zero>(XXMatrix);
   quadrupoleMatrixZero.xy().get<DerivativeOrder::Zero>(XYMatrix);
   quadrupoleMatrixZero.xz().get<DerivativeOrder::Zero>(XZMatrix);
   quadrupoleMatrixZero.yy().get<DerivativeOrder::Zero>(YYMatrix);
   quadrupoleMatrixZero.yz().get<DerivativeOrder::Zero>(YZMatrix);
   quadrupoleMatrixZero.zz().get<DerivativeOrder::Zero>(ZZMatrix);

   // Fill with the first derivative integrals
   quadrupoleMatrixFirst.xx().get<DerivativeOrder::One>(XXderivatives);
   quadrupoleMatrixFirst.xy().get<DerivativeOrder::One>(XYderivatives);
   quadrupoleMatrixFirst.xz().get<DerivativeOrder::One>(XZderivatives);
   quadrupoleMatrixFirst.yy().get<DerivativeOrder::One>(YYderivatives);
   quadrupoleMatrixFirst.yz().get<DerivativeOrder::One>(YZderivatives);
   quadrupoleMatrixFirst.zz().get<DerivativeOrder::One>(ZZderivatives);

   // Get the zeroth and first derivatives of the first element of the first derivative matrix
   quadrupoleMatrixZero.xx().get<DerivativeOrder::Zero>()(0, 0);
   // Get the zeroth and first derivatives of the first element of the first derivative matrix
   quadrupoleMatrixFirst.xx().get<DerivativeOrder::One>()(0, 0).value();
   quadrupoleMatrixFirst.xx().get<DerivativeOrder::One>()(0, 0).derivatives();
 * @endcode
 */
class QuadrupoleMatrix {
 public:
  /**
   * @brief Rule of 6
   * @{
   */
  QuadrupoleMatrix();
  QuadrupoleMatrix(const QuadrupoleMatrix& rhs);
  QuadrupoleMatrix& operator=(const QuadrupoleMatrix& rhs);
  QuadrupoleMatrix(QuadrupoleMatrix&& rhs) noexcept;
  QuadrupoleMatrix& operator=(QuadrupoleMatrix&& rhs) noexcept;
  ~QuadrupoleMatrix();
  //! @}
  /**
   * @brief Resets the quadrupole matrix.
   * @param dimension the number of basis function in the system, i.e. the dimension of the integral matrices.
   */
  void reset(int dimension);
  /**
   * @brief Setters for the xx, xy, xz, yy, yz, and zz components of the matrix.
   * Note: the quadrupole moment integral matrix is symmetric. So, only upper triangle is exposed.
   * @{
   */
  MatrixWithDerivatives& xx() {
    return quadrupoleMatrices_[0];
  }
  MatrixWithDerivatives& xy() {
    return quadrupoleMatrices_[1];
  }
  MatrixWithDerivatives& xz() {
    return quadrupoleMatrices_[2];
  }
  MatrixWithDerivatives& yy() {
    return quadrupoleMatrices_[3];
  }
  MatrixWithDerivatives& yz() {
    return quadrupoleMatrices_[4];
  }
  MatrixWithDerivatives& zz() {
    return quadrupoleMatrices_[5];
  }
  //! @}
  /**
   * @brief Setters for the xx, xy, xz, yy, yz, and zz components of the matrix.
   * Note: the quadrupole moment integral matrix is symmetric. So, only upper triangle is exposed.
   * @{
   */
  const MatrixWithDerivatives& xx() const {
    return quadrupoleMatrices_[0];
  }
  const MatrixWithDerivatives& xy() const {
    return quadrupoleMatrices_[1];
  }
  const MatrixWithDerivatives& xz() const {
    return quadrupoleMatrices_[2];
  }
  const MatrixWithDerivatives& yy() const {
    return quadrupoleMatrices_[3];
  }
  const MatrixWithDerivatives& yz() const {
    return quadrupoleMatrices_[4];
  }
  const MatrixWithDerivatives& zz() const {
    return quadrupoleMatrices_[5];
  }
  //! @}

  /**
   * @brief Returns the components by index. Indices are ordered as xx, xy, xz, yy, yz, zz.
   * @{
   */
  Eigen::MatrixXd& operator[](int index);
  const Eigen::MatrixXd& operator[](int index) const;
  //! @}
  /**
   * @brief Returns the derivative matrices components by index.
   * @{
   */
  MatrixWithDerivatives& at(int index) {
    return quadrupoleMatrices_.at(index);
  }

  const MatrixWithDerivatives& at(int index) const {
    return quadrupoleMatrices_.at(index);
  }
  //! @}
 private:
  std::array<MatrixWithDerivatives, 6> quadrupoleMatrices_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_QUADRUPOLEMATRIX_H
