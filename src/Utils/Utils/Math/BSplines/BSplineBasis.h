/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_BSPLINEBASIS_H
#define BSPLINES_BSPLINEBASIS_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {

namespace BSplines {

/*!
 * Evaluate the B-Spline basis functions that are defined by the knot vector
 */
class BSplineBasis {
 public:
  /*!
   * Evaluates the basis function with the Cox-DeBoor-Mansfield Recurrence Relation
   * - i can range from 0 to n with n being the number of polynomial segements spline
   * - the number of control points or number of basis functions = n+1
   * - n is the number of polynomial segments
   * - the knot vector can be of any degree k
   * - u ranges form 0 to 1
   */
  static double evaluate(int i, int pk, int numberOfSplineSegments, const Eigen::VectorXd& knotVector, double u);
};

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_BSPLINEBASIS_H
