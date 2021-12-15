/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Generator.h"
#include <utility>

namespace Scine {
namespace Utils {

namespace BSplines {

/* Constructor for ControlPolygonGenerator and InterpolationGenerator */
Generator::Generator(const Eigen::MatrixXd& dataPoints, int splineDegree)
  : Generator(dataPoints, static_cast<int>(dataPoints.rows()), splineDegree) {
}

/* Constructor for BSplineFromPenalizedLeastSquaresFit Methods */
Generator::Generator(const Eigen::MatrixXd& dataPoints, int numberOfControlPoints, int splineDegree)
  : dataPoints_(dataPoints),
    p_(splineDegree),
    dim_(static_cast<int>(dataPoints.cols())),
    m_(static_cast<int>(dataPoints.rows()) - 1),
    n_(numberOfControlPoints - 1) {
  assert((p_ >= 1) && "The spline degree p must be at least 1.");
  assert((n_ >= p_) && "The number polynomial segments n must be greater or equal to the spline degree.");
} // namespace BSplines

BSpline Generator::generateBSpline() {
  U_ = generateKnotVector();
  controlPoints_ = generateControlPoints();
  auto spline = BSpline(U_, controlPoints_, p_);
  return spline;
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
