/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ControlPolygonGenerator.h"
#include "GeneratorUtils.h"

namespace Scine {
namespace Utils {

namespace BSplines {

ControlPolygonGenerator::ControlPolygonGenerator(const Eigen::MatrixXd& dataPoints, const unsigned degree,
                                                 const bool uniformKnotVector)
  : Generator(dataPoints, degree), uniformKnotVector_(uniformKnotVector) {
}

Eigen::VectorXd ControlPolygonGenerator::generateKnotVector() {
  if (uniformKnotVector_) {
    return GeneratorUtils::generateKnotVectorByUniformMethod(p_, n_);
  }

  uBar_ = GeneratorUtils::generateParametersByChordLengthMethod(dataPoints_);
  return GeneratorUtils::generateKnotVectorByDeBoorsMethod(p_, n_, uBar_);
}

Eigen::MatrixXd ControlPolygonGenerator::generateControlPoints() {
  return dataPoints_;
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
