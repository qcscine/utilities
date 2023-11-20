/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_FIXEDENDSPENALIZEDLEASTSQUARES_H
#define BSPLINES_FIXEDENDSPENALIZEDLEASTSQUARES_H

#include <Eigen/Core>
#include <Eigen/SVD>

namespace Scine {
namespace Utils {

namespace BSplines {

/*!
 * Generate control points for a B-spline fitted from data points at given coordinates for a known knot vector.
 * \sa FixedEndsPenalizedLeastSquaresGenerator, FixedEndsFixedKnotSplineGenerator
 */
class FixedEndsPenalizedLeastSquares {
 public:
  FixedEndsPenalizedLeastSquares(const Eigen::Ref<const Eigen::MatrixXd>& dataPoints,
                                 const Eigen::VectorXd& dataCoordinates, const Eigen::VectorXd& knotVector,
                                 int splineDegree, bool uniformKnotVector = false, double lambda = 0, int kappa = 2);

  Eigen::MatrixXd calculateControlPoints();

 private:
  /*! eq. (9.63) in the NURBS book */
  Eigen::MatrixXd generateQVectors();

  /*! eq. (9.67) in the NURBS book */
  Eigen::MatrixXd calculateConstantTermsMatrix(const Eigen::MatrixXd& Q);

  /*! eq. (9.66) in the NURBS book */
  Eigen::MatrixXd calculateCoefficientMatrix();
  Eigen::MatrixXd calculateFiniteDifferenceMatrix() const;
  void initializeSolver(const Eigen::MatrixXd& Nmat);
  Eigen::MatrixXd generateControlPointMatrix(const Eigen::MatrixXd& Qmat);

  const Eigen::Ref<const Eigen::MatrixXd>& dataPoints_;
  const Eigen::VectorXd& dataCoordinates_;
  const Eigen::VectorXd& knotVector_;
  int splineDegree_;
  bool uniformKnotVector_;
  double lambda_;
  int kappa_;
  int n_;
  int m_;
  int dim_;

  Eigen::JacobiSVD<Eigen::MatrixXd> svdOfNmat_;
};

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_FIXEDENDSPENALIZEDLEASTSQUARES_H
