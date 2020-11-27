/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_LINEARINTERPOLATOR_H
#define BSPLINES_LINEARINTERPOLATOR_H

#include <Utils/Typenames.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {

namespace BSplines {
class BSpline;

/*!
 * This class generates a bspline interpolating between two given position collections.
 */
class LinearInterpolator {
 public:
  //  static BSpline generate(const Utils::PositionCollection& start, const Utils::PositionCollection& end, int
  //    numberControlPoints);
  static BSpline generate(const Eigen::VectorXd& start, const Eigen::VectorXd& end, int numberControlPoints);
};

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_LINEARINTERPOLATOR_H
