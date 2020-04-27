/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "FixedEndsPenalizedLeastSquares.h"
#include "BSplineBasis.h"
#include "BSplineTools.h"

namespace Scine {
namespace Utils {

namespace BSplines {

FixedEndsPenalizedLeastSquares::FixedEndsPenalizedLeastSquares(const Eigen::Ref<const Eigen::MatrixXd>& dataPoints,
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
  dim_ = static_cast<int>(dataPoints_.cols());
}

Eigen::MatrixXd FixedEndsPenalizedLeastSquares::calculateControlPoints() {
  auto Q = generateQVectors();
  auto Qmat = calculateConstantTermsMatrix(Q);
  auto Nmat = calculateCoefficientMatrix();
  initializeSolver(Nmat);
  return generateControlPointMatrix(Qmat);
}

Eigen::MatrixXd FixedEndsPenalizedLeastSquares::calculateCoefficientMatrix() {
  Eigen::MatrixXd Nmat(m_ - 1, n_ - 1);
  for (int g = 1; g <= m_ - 1; ++g) {
    for (int i = 1; i <= n_ - 1; ++i) {
      Nmat(g - 1, i - 1) = BSplineBasis::evaluate(i, splineDegree_, n_, knotVector_, dataCoordinates_(g));
    }
  }
  return Nmat;
}

Eigen::MatrixXd FixedEndsPenalizedLeastSquares::generateQVectors() {
  Eigen::MatrixXd Q(m_ + 1, dim_);
  Q.setZero();

  for (int g = 1; g <= m_ - 1; ++g) {
    Q.row(g) = dataPoints_.row(g);
    double rec1 = BSplineBasis::evaluate(0, splineDegree_, n_, knotVector_, dataCoordinates_(g));
    double rec2 = BSplineBasis::evaluate(n_, splineDegree_, n_, knotVector_, dataCoordinates_(g));
    for (int l = 0; l < dim_; ++l) {
      Q(g, l) -= rec1 * dataPoints_(0, l);
      Q(g, l) -= rec2 * dataPoints_(m_, l);
    }
  }

  return Q;
}

Eigen::MatrixXd FixedEndsPenalizedLeastSquares::calculateConstantTermsMatrix(const Eigen::MatrixXd& Q) {
  Eigen::MatrixXd Qmat(n_ - 1, dim_);
  for (int i = 1; i <= n_ - 1; ++i) {
    Eigen::VectorXd accum(dim_, 1);
    accum.setZero();

    for (int g = 1; g <= m_ - 1; ++g) {
      accum += BSplineBasis::evaluate(i, splineDegree_, n_, knotVector_, dataCoordinates_(g)) * Q.row(g);
    }
    Qmat.row(i - 1) = accum;
  }
  return Qmat;
}

Eigen::MatrixXd FixedEndsPenalizedLeastSquares::calculateFiniteDifferenceMatrix() {
  assert(kappa_ <= n_);
  assert(uniformKnotVector_ && "For penalization, the knot vector must be equidistant according to Eilers. The "
                               "equations for non-equidistant knot vectors is not yet implemented.");

  Eigen::MatrixXd DeltaMat(n_ - 1 - kappa_, n_ - 1);
  DeltaMat.setZero();

  for (int i = 0; i <= n_ - kappa_ - 2; ++i) {
    for (int j = 0; j <= n_ - 2; ++j) {
      DeltaMat(i, j) = BSplineTools::differenceOperator(i, j, kappa_);
    }
  }

  return DeltaMat;
}

void FixedEndsPenalizedLeastSquares::initializeSolver(const Eigen::MatrixXd& Nmat) {
  Eigen::MatrixXd m = Nmat.transpose() * Nmat;

  if (lambda_ != 0) {
    auto deltaMat = calculateFiniteDifferenceMatrix();
    m += lambda_ * (deltaMat.transpose() * deltaMat);
  }

  svdOfNmat_.compute(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
}

Eigen::MatrixXd FixedEndsPenalizedLeastSquares::generateControlPointMatrix(const Eigen::MatrixXd& Qmat) {
  Eigen::MatrixXd cp(n_ + 1, dim_);

  Eigen::MatrixXd Pinternal(n_ - 1, dim_);
  Pinternal = svdOfNmat_.solve(Qmat);

  // construct control point vector P_
  cp.row(0) = dataPoints_.row(0);
  for (int i = 1; i <= n_ - 1; ++i) {
    cp.row(i) = Pinternal.row(i - 1);
  }
  cp.row(n_) = dataPoints_.row(m_);
  return cp;
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
