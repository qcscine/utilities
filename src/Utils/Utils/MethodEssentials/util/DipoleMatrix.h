/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_DIPOLEMATRIX_H
#define UTILS_DIPOLEMATRIX_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {

class InvalidDimensionSizeException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "Invalid dimension index. Possible dimensions: 0 = x, 1 = y, 2 = z.";
  }
};

/**
 * @class DipoleMatrix @file DipoleMatrix.h
 * This class stores the dipole integrals <mu|r|nu> in their three components x, y, z.
 */
class DipoleMatrix {
 public:
  /**
   * @brief Resets the dipole matrix.
   * @param dimension the number of basis function in the system, i.e. the dimension of the integral matrices.
   */
  void reset(int dimension);
  /**
   * @brief Getters and setters for the x, y and z components of the matrix.
   * @{
   */
  Eigen::MatrixXd& x();
  Eigen::MatrixXd x() const;
  Eigen::MatrixXd& y();
  Eigen::MatrixXd y() const;
  Eigen::MatrixXd& z();
  Eigen::MatrixXd z() const;
  /**
   * @}
   */
  /**
   * @brief Returns the components by index. 0 is x, 1 is y and 2 is z.
   */
  Eigen::MatrixXd& operator[](int index);
  Eigen::MatrixXd operator[](int index) const;
  /**
   * @brief Updates entries after displacement of one atom.
   * This method updates all the integrals needed after the displacement of one atom.
   * This is particularly convenient for the semi-numerical calculation of the dipole gradient.
   */
  // TODO: Implement the block update.
 private:
  Eigen::MatrixXd xDipoleMatrix_;
  Eigen::MatrixXd yDipoleMatrix_;
  Eigen::MatrixXd zDipoleMatrix_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_DIPOLEMATRIX_H
