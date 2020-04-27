/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LinearInterpolator.h"
#include "BSpline.h"
#include "InterpolationGenerator.h"

namespace Scine {
namespace Utils {

namespace BSplines {

// BSpline LinearInterpolator::generate(const Utils::PositionCollection& start, const Utils::PositionCollection& end,
//                                     int numberControlPoints) {
//  Eigen::VectorXd startVector = Eigen::Map<const Eigen::VectorXd>(start.data(), start.cols()*start.rows());
//  Eigen::VectorXd endVector = Eigen::Map<const Eigen::VectorXd>(end.data(), end.cols()*end.rows());
//  return generate(startVector, endVector, numberControlPoints);
//}

BSpline LinearInterpolator::generate(const Eigen::VectorXd& start, const Eigen::VectorXd& end, int numberControlPoints) {
  assert(numberControlPoints >= 2);
  assert(start.size() == end.size());

  int numberInterpolationPoints = numberControlPoints;
  double delta = 1.0 / (numberInterpolationPoints - 1);
  Eigen::MatrixXd interpolationPoints(numberInterpolationPoints, start.size());
  for (int i = 0; i < numberInterpolationPoints; ++i) {
    double endFraction = i * delta;
    interpolationPoints.row(i) = (1.0 - endFraction) * start + endFraction * end;
  }
  InterpolationGenerator generator(interpolationPoints);
  return generator.generateBSpline();
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
