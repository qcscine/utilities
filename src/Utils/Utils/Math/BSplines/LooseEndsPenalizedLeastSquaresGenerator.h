/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_LOOSEENDSPENALIZEDLEASTSQUARESGENERATOR_H
#define BSPLINES_LOOSEENDSPENALIZEDLEASTSQUARESGENERATOR_H

#include "BSpline.h"
#include "Generator.h"
#include <Eigen/SVD>

namespace Scine {
namespace Utils {

namespace BSplines {

/*!
 * Generates a B-Spline curve with a specified number of control points approximating the dat a points with loose ends.
 * „Loose ends" means that the first and last control points do not need to coincide with the first and last data point
 * of the set. Thus, all n+1 control points are optimized. Still, the resulting B-Spline curve passes through the first
 * and last control point since the first and last basis functions are equal to unity at the ends of the domain
 * (N(u=0)=1, N(u=1)=1, therefore C(u=0)=P0 and C(u=1)=Plast). Penalization can be turned on by specifying a lambda
 * value > 0 and difference orders can be specified by kappa.
 */
class LooseEndsPenalizedLeastSquaresGenerator : public Generator {
 public:
  LooseEndsPenalizedLeastSquaresGenerator(const Eigen::Ref<const Eigen::MatrixXd>& dataPoints,
                                          int numberOfControlPoints, int splineDegree = 3,
                                          bool uniformKnotVector = false, double lambda = 0, int kappa = 2);

 private:
  Eigen::VectorXd generateKnotVector() override;
  Eigen::MatrixXd generateControlPoints() override;

  bool uniformKnotVector_;
  double lambda_;
  int kappa_;
};

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_LOOSEENDSPENALIZEDLEASTSQUARESGENERATOR_H
