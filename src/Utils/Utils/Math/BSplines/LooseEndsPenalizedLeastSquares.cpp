/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LooseEndsPenalizedLeastSquares.h"
#include "BSplineBasis.h"
#include "BSplineTools.h"

namespace Scine {
namespace Utils {

namespace BSplines {

LooseEndsPenalizedLeastSquares::LooseEndsPenalizedLeastSquares(const Eigen::Ref<const Eigen::MatrixXd>& dataPoints,
                                                               const Eigen::VectorXd& dataCoordinates,
                                                               const Eigen::VectorXd& knotVector, int splineDegree,
                                                               bool uniformKnotVector, double lambda, int kappa)
  : dataPoints_(dataPoints),
    dataCoordinates_(dataCoordinates),
    knotVector_(knotVector),
    splineDegree_(splineDegree),
    uniformKnotVector_(uniformKnotVector),
    lambda_(lambda),
    kappa_(kappa) {
  assert(dataPoints.rows() == dataCoordinates.size());
  n_ = static_cast<int>(knotVector.size()) - splineDegree_ - 2;
  m_ = static_cast<int>(dataPoints_.rows()) - 1;
}

Eigen::MatrixXd LooseEndsPenalizedLeastSquares::calculateControlPoints() {
  // generateQVectors();
  // calculateConstantTermsMatrix();
  auto Nmat = calculateCoefficientMatrix();
  initializeSolver(Nmat);
  return generateControlPointMatrix(Nmat);
}

Eigen::MatrixXd LooseEndsPenalizedLeastSquares::calculateCoefficientMatrix() {
  Eigen::MatrixXd Nmat(m_ + 1, n_ + 1);
  for (int g = 0; g <= m_; ++g) {
    for (int i = 0; i <= n_; ++i) {
      Nmat(g, i) = BSplineBasis::evaluate(i, splineDegree_, n_, knotVector_, dataCoordinates_(g));
    }
  }
  return Nmat;
}

Eigen::MatrixXd LooseEndsPenalizedLeastSquares::calculateFiniteDifferenceMatrix() const {
  assert(kappa_ <= n_);
  assert(uniformKnotVector_ && "For penalization, the knot vector must be equidistant according to Eilers. The "
                               "equations for non-equidistant knot vectors is not yet implemented.");

  Eigen::MatrixXd DeltaMat_(n_ + 1 - kappa_, n_ + 1);
  DeltaMat_.setZero();

  for (int i = 0; i <= n_ - kappa_; ++i) {
    for (int j = 0; j <= n_; ++j) {
      DeltaMat_(i, j) = BSplineTools::differenceOperator(i, j, kappa_);
    }
  }

  return DeltaMat_;
}

void LooseEndsPenalizedLeastSquares::initializeSolver(const Eigen::MatrixXd& Nmat) {
  Eigen::MatrixXd m = Nmat.transpose() * Nmat;

  if (lambda_ != 0) {
    auto deltaMat = calculateFiniteDifferenceMatrix();
    m += lambda_ * (deltaMat.transpose() * deltaMat);
  }

  svdOfNmat_.compute(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
}

Eigen::MatrixXd LooseEndsPenalizedLeastSquares::generateControlPointMatrix(const Eigen::MatrixXd& Nmat) {
  return svdOfNmat_.solve(Nmat.transpose() * dataPoints_);
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
