/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_CONTROLPOLYGONGENERATOR_H
#define BSPLINES_CONTROLPOLYGONGENERATOR_H

#include "Generator.h"

namespace Scine {
namespace Utils {

namespace BSplines {

/*! Generates a B-spline curve from a set of data points that are used as control points.
 * */
class ControlPolygonGenerator : public Generator {
 public:
  explicit ControlPolygonGenerator(const Eigen::MatrixXd& dataPoints, unsigned degree = 3, bool uniformKnotVector = false);

 private:
  Eigen::VectorXd generateKnotVector() override;
  Eigen::MatrixXd generateControlPoints() override;

  bool uniformKnotVector_;
};

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_CONTROLPOLYGONGENERATOR_H
