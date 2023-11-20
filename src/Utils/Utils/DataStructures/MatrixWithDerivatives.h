/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATRIXDER_H
#define UTILS_MATRIXDER_H

#include <Utils/Math/AutomaticDifferentiation/TypeDefinitions.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {

// clang-format off

/*!
 * @class MatrixWithDerivatives @file MatrixWithDerivatives.h
 * @brief Container for a matrix and the derivatives of its elements.
 * This class contains 3 matrices:
 * val: contains the 0th order derivative of the elements.
 * der: contains the 1th order derivative of the elements.
 *      It is thus a Eigen::Matrix of AutomaticDifferentiation::First3D elements.
 * hes: contains the 2th order derivative of the elements.
 *      It is thus a Eigen::Matrix of AutomaticDifferentiation::Second3D elements.
 */
class MatrixWithDerivatives{
public:
  using Der0 = double;
  using Der1 = AutomaticDifferentiation::First3D;
  using Der2 = AutomaticDifferentiation::Second3D;
  using Matrix0 = Eigen::MatrixXd;
  using Matrix1 = Eigen::Matrix<Der1, Eigen::Dynamic, Eigen::Dynamic>;
  using Matrix2 = Eigen::Matrix<Der2, Eigen::Dynamic, Eigen::Dynamic>;
  template< DerivativeOrder O > struct MatrixType { using MType = Matrix0; }; //Default
  template<DerivativeOrder o> using Matrix = typename MatrixType<o>::MType;

  MatrixWithDerivatives& operator+(const MatrixWithDerivatives& rhs);
  MatrixWithDerivatives& operator+=(const MatrixWithDerivatives& rhs);
  MatrixWithDerivatives& operator-(const MatrixWithDerivatives& rhs);
  MatrixWithDerivatives& operator-=(const MatrixWithDerivatives& rhs);

  explicit MatrixWithDerivatives(int rows=0, int cols=0) {
    order_ = DerivativeOrder::Zero;
    setDimension(rows, cols);
  }

  void setOrder(DerivativeOrder o) {
    order_ = o;
  }

  //TODO: ASSERTS IF WRONG ORDER IS CALLED
  Der0& v0(int i1, int i2) { return val(i1, i2); }
  Der1& v1(int i1, int i2) { return der(i1, i2); }
  Der2& v2(int i1, int i2) { return hes(i1, i2); }
  const Der0& v0(int i1, int i2) const { return val(i1, i2); }
  const Der1& v1(int i1, int i2) const { return der(i1, i2); }
  const Der2& v2(int i1, int i2) const { return hes(i1, i2); }
  template<DerivativeOrder O> Matrix<O>& get();
  template<DerivativeOrder O> const Matrix<O>& get() const;
  Matrix0& get0() { return val; }
  Matrix1& get1() { return der; }
  Matrix2& get2() { return hes; }
  const Matrix0& get0() const { return val; }
  const Matrix1& get1() const { return der; }
  const Matrix2& get2() const { return hes; }

  /*! Get a copy of the underlying matrix, without derivatives. */
  Eigen::MatrixXd getMatrixXd() const;

  /*! Alias for setDimension for the case where rows = cols */
  void reset(int dimension) {
    setOrder(DerivativeOrder::Zero);
    setDimension(dimension, dimension);
  }
  /*! Initializes the members with dimensions given as parameters. */
  void setDimension(int rows, int cols);

  /*! Sets the base matrix to m.
      Implicitly calls setDimension with correct dimensions. */
  void setBaseMatrix(const Eigen::MatrixXd& m);

  double getValue(int i, int j) const {
    if (order_ == DerivativeOrder::Zero) {
      return v0(i, j);
}
    if (order_ == DerivativeOrder::One) {
      return v1(i, j).value();
}
    return v2(i, j).value();
  }

private:
  DerivativeOrder order_;
  int nCols_ = 0, nRows_ = 0;
  Matrix0 val;
  Matrix1 der;
  Matrix2 hes;
};

template<> struct MatrixWithDerivatives::MatrixType<DerivativeOrder::One> { using MType = Matrix1; };
template<> struct MatrixWithDerivatives::MatrixType<DerivativeOrder::Two> { using MType = Matrix2; };
template<> inline MatrixWithDerivatives::Matrix0& MatrixWithDerivatives::get<DerivativeOrder::Zero>() { return get0(); }
template<> inline MatrixWithDerivatives::Matrix1& MatrixWithDerivatives::get<DerivativeOrder::One>() { return get1(); }
template<> inline MatrixWithDerivatives::Matrix2& MatrixWithDerivatives::get<DerivativeOrder::Two>() { return get2(); }
template<> inline const MatrixWithDerivatives::Matrix0& MatrixWithDerivatives::get<DerivativeOrder::Zero>() const { return get0(); }
template<> inline const MatrixWithDerivatives::Matrix1& MatrixWithDerivatives::get<DerivativeOrder::One>() const { return get1(); }
template<> inline const MatrixWithDerivatives::Matrix2& MatrixWithDerivatives::get<DerivativeOrder::Two>() const { return get2(); }

// clang-format on

} // namespace Utils
} // namespace Scine
#endif // UTILS_MATRIXDER_H
