/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_OCTUPOLEMATRIX_H
#define UTILSOS_OCTUPOLEMATRIX_H

#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <array>

namespace Scine {
namespace Utils {

/**
 * @class OctupoleMatrix @file OctupoleMatrix.h
 * @brief Class representing octupole matrix integrals with derivatives up to 2nd order.
 *
 * This class stores the octupole integrals \f$\langle mu|o_{rr}|nu \rangle\f$
 * The class stores only the upper diagonal of the octupole moment tensor:
 * xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz.
 * The derivatives are accessible through the .derivatives() method of the
 * underlying AutomaticDifferentiation::First3D/Second3D class.
 *
 * @code{cpp}
   int arbitraryDimension = 10;
   OctupoleMatrix octupoleMatrixZero, octupoleMatrixFirst;

   // Fill OctupoleMatrix with random values and derivatives.
    Eigen::MatrixXd XXXMatrix, YYYMatrix, ZZZMatrix, XXYMatrix, XXZMatrix, XYZMatrix,
                    XYYMatrix, XZZMatrix, YYZMatrix, YZZMatrix;
   Eigen::Matrix<AutomaticDifferentiation::First3D, Eigen::Dynamic, Eigen::Dynamic> XXXderivatives,
     YYYderivatives, ZZZderivatives, XXYderivatives, XXZderivatives, XYZderivatives,
     XYYderivatives, XZZderivatives, YYZderivatives, YZZderivatives;

   XXXMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   YYYMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   ZZZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   XXYMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   XXZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   XYZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   XYYMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   XZZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   YYZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   YZZMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
   XYYderivatives.resize(arbitraryDimension, arbitraryDimension);
   XZZderivatives.resize(arbitraryDimension, arbitraryDimension);
   YYZderivatives.resize(arbitraryDimension, arbitraryDimension);
   YZZderivatives.resize(arbitraryDimension, arbitraryDimension);
   for (int i = 0; i < arbitraryDimension; ++i) {
     for (int j = 0; j < arbitraryDimension; ++j) {
       XXXderivatives(i, j) = {XXXMatrix(i, j), Eigen::Vector3d::Random()};
       YYYderivatives(i, j) = {YYYMatrix(i, j), Eigen::Vector3d::Random()};
       ZZZderivatives(i, j) = {ZZZMatrix(i, j), Eigen::Vector3d::Random()};
       XXYderivatives(i, j) = {XXYMatrix(i, j), Eigen::Vector3d::Random()};
       XXZderivatives(i, j) = {XXZMatrix(i, j), Eigen::Vector3d::Random()};
       XYZderivatives(i, j) = {XYZMatrix(i, j), Eigen::Vector3d::Random()};
       XYYderivatives(i, j) = {XYYMatrix(i, j), Eigen::Vector3d::Random()};
       XZZderivatives(i, j) = {XZZMatrix(i, j), Eigen::Vector3d::Random()};
       YYZderivatives(i, j) = {YYZMatrix(i, j), Eigen::Vector3d::Random()};
       YZZderivatives(i, j) = {YZZMatrix(i, j), Eigen::Vector3d::Random()};
     }
   }
   // Fill with the zero-th derivative integrals
   octupoleMatrixZero.xxx().get<derivOrder::zero>(XXXXMatrix);
   octupoleMatrixZero.yyy().get<derivOrder::zero>(YYYMatrix);
   octupoleMatrixZero.zzz().get<derivOrder::zero>(ZZZMatrix);
   octupoleMatrixZero.xxy().get<derivOrder::zero>(XXYMatrix);
   octupoleMatrixZero.xxz().get<derivOrder::zero>(XXZMatrix);
   octupoleMatrixZero.xyz().get<derivOrder::zero>(XYZMatrix);
   octupoleMatrixZero.xyy().get<derivOrder::zero>(XYYMatrix);
   octupoleMatrixZero.xzz().get<derivOrder::zero>(XZZMatrix);
   octupoleMatrixZero.yyz().get<derivOrder::zero>(YYZMatrix);
   octupoleMatrixZero.yzz().get<derivOrder::zero>(YZZMatrix);

   // Fill with the first derivative integrals
   octupoleMatrixFirst.xxx().get<derivOrder::one>(XXXderivatives);
   octupoleMatrixFirst.yyy().get<derivOrder::one>(YYYderivatives);
   octupoleMatrixFirst.zzz().get<derivOrder::one>(ZZZderivatives);
   octupoleMatrixFirst.xxy().get<derivOrder::one>(XXYderivatives);
   octupoleMatrixFirst.xxz().get<derivOrder::one>(XXZderivatives);
   octupoleMatrixFirst.xyz().get<derivOrder::one>(XYZderivatives);
   octupoleMatrixFirst.xyy().get<derivOrder::one>(XYYderivatives);
   octupoleMatrixFirst.xzz().get<derivOrder::one>(XZZderivatives);
   octupoleMatrixFirst.yyz().get<derivOrder::one>(YYZderivatives);
   octupoleMatrixFirst.yzz().get<derivOrder::one>(YZZderivatives);

   // Get the zeroth and first derivatives of the first element of the first derivative matrix
   octupoleMatrixZero.xxx().get<derivOrder::zero>()(0, 0);
   // Get the zeroth and first derivatives of the first element of the first derivative matrix
   octupoleMatrixFirst.xxx().get<derivOrder::one>()(0, 0).value();
   // or
   AutomaticDifferentiation::getValue3DAsDouble(octupoleMatrixFirst.xxx().get<derivOrder::one>()(0, 0));
   octupoleMatrixFirst.xxx().get<derivOrder::one>()(0, 0).derivatives();
 * @endcode
 */
class OctupoleMatrix {
 public:
  /**
   * @brief Rule of 6
   * @{
   */
  OctupoleMatrix();
  OctupoleMatrix(const OctupoleMatrix& rhs);
  OctupoleMatrix& operator=(const OctupoleMatrix& rhs);
  OctupoleMatrix(OctupoleMatrix&& rhs) noexcept;
  OctupoleMatrix& operator=(OctupoleMatrix&& rhs) noexcept;
  ~OctupoleMatrix();
  //! @}
  /**
   * @brief Resets the quadrupole matrix.
   * @param dimension the number of basis function in the system, i.e. the dimension of the integral matrices.
   */
  void reset(int dimension);
  /**
   * @brief Setters for the xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz and zzz components of the tensor.
   * Note: the quadrupole moment integral matrix is symmetric. So, only upper triangle is exposed.
   * @{
   */
  MatrixWithDerivatives& xxx() {
    return octupoleMatrices_[0];
  }
  MatrixWithDerivatives& xxy() {
    return octupoleMatrices_[1];
  }
  MatrixWithDerivatives& xxz() {
    return octupoleMatrices_[2];
  }
  MatrixWithDerivatives& xyy() {
    return octupoleMatrices_[3];
  }
  MatrixWithDerivatives& xyz() {
    return octupoleMatrices_[4];
  }
  MatrixWithDerivatives& xzz() {
    return octupoleMatrices_[5];
  }
  MatrixWithDerivatives& yyy() {
    return octupoleMatrices_[6];
  }
  MatrixWithDerivatives& yyz() {
    return octupoleMatrices_[7];
  }
  MatrixWithDerivatives& yzz() {
    return octupoleMatrices_[8];
  }
  MatrixWithDerivatives& zzz() {
    return octupoleMatrices_[9];
  }
  //! @}
  /**
   * @brief Getters for the xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz and zzz components of the matrix.
   * Note: the quadrupole moment integral matrix is symmetric. So, only upper triangle is exposed.
   * @{
   */
  const MatrixWithDerivatives& xxx() const {
    return octupoleMatrices_[0];
  }
  const MatrixWithDerivatives& xxy() const {
    return octupoleMatrices_[1];
  }
  const MatrixWithDerivatives& xxz() const {
    return octupoleMatrices_[2];
  }
  const MatrixWithDerivatives& xyy() const {
    return octupoleMatrices_[3];
  }
  const MatrixWithDerivatives& xyz() const {
    return octupoleMatrices_[4];
  }
  const MatrixWithDerivatives& xzz() const {
    return octupoleMatrices_[5];
  }
  const MatrixWithDerivatives& yyy() const {
    return octupoleMatrices_[6];
  }
  const MatrixWithDerivatives& yyz() const {
    return octupoleMatrices_[7];
  }
  const MatrixWithDerivatives& yzz() const {
    return octupoleMatrices_[8];
  }
  const MatrixWithDerivatives& zzz() const {
    return octupoleMatrices_[9];
  }
  //! @}

  /**
   * @brief Returns the components by index. Indices are ordered as xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz.
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
    return octupoleMatrices_.at(index);
  }

  const MatrixWithDerivatives& at(int index) const {
    return octupoleMatrices_.at(index);
  }
  //! @}
 private:
  std::array<MatrixWithDerivatives, 10> octupoleMatrices_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_OCTUPOLEMATRIX_H
