/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_FIXEDENDSPENALIZEDSQUARESGENERATOR_H
#define BSPLINES_FIXEDENDSPENALIZEDSQUARESGENERATOR_H

#include "Generator.h"

namespace Scine {
namespace Utils {

namespace BSplines {

/*!
 * Generates a B-Spline curve with a specified number of control points approximating the data points with fixed ends.
 * „Fixed ends" means that the first and last control points are constrained to coincide with the first and last data
 * point of the set. (see the NURBS book by Piegl 1997) Thus, all but the first and last control points (n-1 in total)
 * are optimized. Still, the resulting B-Splines curve passes through the first and last control point since the first
 * and last basis functions are equal to unity at the ends of the domain (N(u=0)=1, N(u=1)=1, therefore C(u=0)=P_0=R_0
 * and C(u=1)=P_last=R_last). Penalization can be turned on by specifying a lambda value > 0 and difference orders can
 * be specified by kappa. (see the NURBS book by Piegl 1997)
 */
class FixedEndsPenalizedLeastSquaresGenerator : public Generator {
 public:
  FixedEndsPenalizedLeastSquaresGenerator(const Eigen::Ref<const Eigen::MatrixXd>& dataPoints,
                                          int numberOfControlPoints, int splineDegree = 3,
                                          bool uniformKnotVector = false, double lambda = 0, int kappa = 2);

 private:
  Eigen::VectorXd generateKnotVector() override;
  Eigen::MatrixXd generateControlPoints() override;

  bool uniformKnotVector_;
  double lambda_;
  int kappa_;
  Eigen::MatrixXd Q_, Qmat_, Nmat_, DeltaMat_;
};

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_FIXEDENDSPENALIZEDSQUARESGENERATOR_H
