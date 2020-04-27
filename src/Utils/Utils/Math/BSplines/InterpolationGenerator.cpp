/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "InterpolationGenerator.h"
#include "BSplineBasis.h"
#include "GeneratorUtils.h"

namespace Scine {
namespace Utils {

namespace BSplines {

InterpolationGenerator::InterpolationGenerator(const Eigen::Ref<const Eigen::MatrixXd>& dataPoints,
                                               const unsigned degree, const bool uniformKnotVector)
  : Generator(dataPoints, degree), uniformKnotVector_(uniformKnotVector) {
}

Eigen::VectorXd InterpolationGenerator::generateKnotVector() {
  if (uniformKnotVector_) {
    return GeneratorUtils::generateKnotVectorByUniformMethod(p_, n_);
  }
  else {
    // generateParametersByChordLengthMethod();
    uBar_ = GeneratorUtils::generateParametersByCentripetalMethod(dataPoints_);
    return GeneratorUtils::generateKnotVectorByKnotAveraging(p_, n_, uBar_);
  }
}

Eigen::MatrixXd InterpolationGenerator::generateControlPoints() {
  calculateCoefficientMatrix();
  initializeSolver();
  return generateControlPointMatrix();
}

void InterpolationGenerator::calculateCoefficientMatrix() {
  Nmat_.resize(n_ + 1, n_ + 1);
  for (int g = 0; g <= m_; ++g) {
    for (int i = 0; i <= n_; ++i) {
      Nmat_(g, i) = BSplineBasis::evaluate(i, p_, n_, U_, uBar_(g));
    }
  }
}

void InterpolationGenerator::initializeSolver() {
  // svdOfNmat_.compute(Nmat_);
  svdOfNmat_.compute(Nmat_, Eigen::ComputeThinU | Eigen::ComputeThinV);
}

Eigen::MatrixXd InterpolationGenerator::generateControlPointMatrix() {
  return svdOfNmat_.solve(dataPoints_);
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
