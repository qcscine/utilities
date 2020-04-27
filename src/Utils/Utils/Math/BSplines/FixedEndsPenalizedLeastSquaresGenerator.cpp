/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "FixedEndsPenalizedLeastSquaresGenerator.h"
#include "BSplineBasis.h"
#include "FixedEndsPenalizedLeastSquares.h"
#include "GeneratorUtils.h"

namespace Scine {
namespace Utils {

namespace BSplines {

FixedEndsPenalizedLeastSquaresGenerator::FixedEndsPenalizedLeastSquaresGenerator(const Eigen::Ref<const Eigen::MatrixXd>& dataPoints,
                                                                                 int numberOfControlPoints,
                                                                                 int splineDegree, bool uniformKnotVector,
                                                                                 double lambda, int kappa)
  : Generator(dataPoints, numberOfControlPoints, splineDegree),
    uniformKnotVector_(uniformKnotVector),
    lambda_(lambda),
    kappa_(kappa) {
  assert((m_ > n_) && "The number of data segments m has to be greater than "
                      "the number of requested polynomial segments n.");
}

Eigen::VectorXd FixedEndsPenalizedLeastSquaresGenerator::generateKnotVector() {
  if (!uniformKnotVector_) {
    // generateParametersByChordLengthMethod();
    uBar_ = GeneratorUtils::generateParametersByCentripetalMethod(dataPoints_);

    // generateKnotVectorByKnotAveraging();
    return GeneratorUtils::generateKnotVectorByDeBoorsMethod(p_, n_, uBar_);
  }
  else {
    uBar_ = GeneratorUtils::generateParametersByEquallySpacedMethod(dataPoints_);
    return GeneratorUtils::generateKnotVectorByUniformMethod(p_, n_);
  }
}

Eigen::MatrixXd FixedEndsPenalizedLeastSquaresGenerator::generateControlPoints() {
  FixedEndsPenalizedLeastSquares ls(dataPoints_, uBar_, U_, p_, uniformKnotVector_, lambda_, kappa_);
  return ls.calculateControlPoints();
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
