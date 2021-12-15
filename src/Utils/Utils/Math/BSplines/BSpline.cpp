/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "BSpline.h"
#include "BSplineBasis.h"
#include "Coefficients.h"
#include "Exceptions.h"
#include <iostream>

namespace Scine {
namespace Utils {

namespace BSplines {

BSpline::BSpline() {
  // Give the knot vector the minimal number of elements it needs.
  Eigen::VectorXd knotVector(2);
  knotVector << 0, 1;
  Uk_.push_back(std::move(knotVector));
  // Create a minimal control point matrix
  Eigen::MatrixXd controlPoints(1, 0);
  Pk_.push_back(std::move(controlPoints));
}

BSpline::BSpline(Eigen::VectorXd knotVector, Eigen::MatrixXd controlPoints, int degree)
  : p_(degree), n_(static_cast<int>(controlPoints.rows()) - 1), dim_(static_cast<int>(controlPoints.cols())) {
  assert(p_ >= 1 && "A B-spline curve must have a degree of at least 1");
  assert(n_ >= p_ && "A B-spline curve must have more control points than its own degree.");

  // Reserve the maximal size, to avoid possible dangling references when the derivatives must be calculated
  Uk_.reserve(p_ + 1);
  Pk_.reserve(p_ + 1);

  Uk_.push_back(std::move(knotVector));
  Pk_.push_back(std::move(controlPoints));
}

BSpline BSpline::getDerivativeBSpline(int derivativeOrder) {
  assert(0 <= derivativeOrder && derivativeOrder < p_);
  return BSpline(deriveAndGetKnotVector(derivativeOrder), deriveAndGetControlPointMatrix(derivativeOrder), p_ - derivativeOrder);
}

void BSpline::reverse() {
  for (auto& knotVector : Uk_) {
    Eigen::VectorXd tempVec(knotVector);
    knotVector = 1.0 - tempVec.reverse().array();
  }

  for (auto& controlPointMatrix : Pk_) {
    controlPointMatrix.colwise().reverseInPlace();
  }
}

BSpline BSpline::reversed() const {
  auto bspline = *this;
  bspline.reverse();
  return bspline;
}

bool BSpline::isClampedAndNormed() const {
  Eigen::VectorXd head;
  Eigen::VectorXd tail;
  head = Eigen::VectorXd::Zero(p_ + 1);
  tail = Eigen::VectorXd::Ones(p_ + 1);

  return Uk_[0].head(p_ + 1).isApprox(head) && Uk_[0].tail(p_ + 1).isApprox(tail);
}

const Eigen::VectorXd& BSpline::getKnotVector() const {
  return Uk_[0];
}

const Eigen::VectorXd& BSpline::getKnotVectorDerivative(int derivativeOrder) const {
  if (derivativeOrder > highestCalculatedDerivative_) {
    calculateDerivatives(derivativeOrder);
  }

  return Uk_[derivativeOrder];
}

const Eigen::VectorXd& BSpline::deriveAndGetKnotVector(int derivativeOrder) {
  assert((0 <= derivativeOrder && derivativeOrder < p_) && "This derivative does not exist.");
  if (derivativeOrder > highestCalculatedDerivative_) {
    calculateDerivatives(derivativeOrder);
  }
  return getKnotVectorDerivative(derivativeOrder);
}

const Eigen::MatrixXd& BSpline::getControlPointMatrix() const {
  return getControlPointMatrix(0);
}

const Eigen::MatrixXd& BSpline::getControlPointMatrix(int k) const {
  if (k > highestCalculatedDerivative_) {
    calculateDerivatives(k);
  }

  return Pk_[k];
}

const Eigen::MatrixXd& BSpline::deriveAndGetControlPointMatrix(int derivativeOrder) {
  assert((0 <= derivativeOrder && derivativeOrder < p_) && "This derivative does not exist.");
  if (derivativeOrder > highestCalculatedDerivative_) {
    calculateDerivatives(derivativeOrder);
  }
  return getControlPointMatrix(derivativeOrder);
}

void BSpline::deriveKnotVectors(int highestDerivativeToCalculate) const {
  assert((highestDerivativeToCalculate >= 0) && (highestDerivativeToCalculate <= p_));

  for (auto k = static_cast<int>(Uk_.size()); k <= highestDerivativeToCalculate; ++k) {
    /* make the dth derivative knot vector from the previous d-1th derivative knot vector segment by dropping the
     * first and last segment */
    Eigen::VectorXd segmentOfPreviousOrderKnotVector(Uk_[k - 1].segment(1, Uk_[k - 1].size() - 2));
    Uk_.push_back(segmentOfPreviousOrderKnotVector);
  }
}

void BSpline::deriveControlPointMatrices(int highestDerivativeToCalculate) const {
  assert((highestDerivativeToCalculate >= 0) && (highestDerivativeToCalculate <= p_));

  for (auto k = static_cast<int>(Pk_.size()); k <= highestDerivativeToCalculate; ++k) {
    Eigen::MatrixXd Pnew(n_ - k + 1, dim_);
    for (int i = 0; i <= n_ - k; ++i) {
      Pnew.row(i) = deriveControlPoint(i, k);
    }
    Pk_.push_back(Pnew);
  }
}

void BSpline::calculateDerivatives(int highestDerivativeToCalculate) const {
  assert(highestDerivativeToCalculate >= 0);
  if (highestDerivativeToCalculate > p_) {
    throw std::runtime_error("Invalid derivative order");
  }

  deriveKnotVectors(highestDerivativeToCalculate);
  deriveControlPointMatrices(highestDerivativeToCalculate);
  highestCalculatedDerivative_ = highestDerivativeToCalculate;
}

Eigen::VectorXd BSpline::deriveControlPoint(int i, int highestDerivativeToCalculate) const {
  if (highestDerivativeToCalculate == 0) {
    return Pk_[0].row(i);
  }

  if (Uk_[0](i + p_ + 1) == Uk_[0](i + highestDerivativeToCalculate)) {
    return Eigen::VectorXd(Eigen::VectorXd::Zero(dim_));
  }

  return double(p_ - highestDerivativeToCalculate + 1) / (Uk_[0](i + p_ + 1) - Uk_[0](i + highestDerivativeToCalculate)) *
         (deriveControlPoint(i + 1, highestDerivativeToCalculate - 1) -
          deriveControlPoint(i, highestDerivativeToCalculate - 1));
}

Eigen::VectorXd BSpline::evaluateNaive(double u, int derivativeOrder) const {
  assert((0.0 <= u) && (u <= 1.0));

  if (derivativeOrder > highestCalculatedDerivative_) {
    calculateDerivatives(derivativeOrder);
  }

  Eigen::VectorXd C(dim_);
  C.setZero();
  for (int i = 0; i <= n_ - derivativeOrder; ++i) {
    C += BSplineBasis::evaluate(i, p_ - derivativeOrder, n_ - derivativeOrder, Uk_[derivativeOrder], u) *
         Pk_[derivativeOrder].row(i);
  }
  return C;
}

Eigen::VectorXd BSpline::operator()(double u) const {
  return evaluate(u, 0);
}

Eigen::VectorXd BSpline::calculateBSplineCoefficientVector(double u, int derivativeOrder) const {
  return calculateBSplineCoefficients(u, derivativeOrder).fullCoefficientVector();
}

Coefficients BSpline::calculateBSplineCoefficients(double u, int derivativeOrder) const {
  if (derivativeOrder > highestCalculatedDerivative_) {
    calculateDerivatives(derivativeOrder);
  }

  int l = findIdxOfLowerOrEqualDomainKnot(u, derivativeOrder);
  Eigen::VectorXd allCoefficients = deBoorCoefficient(u, l, p_ - derivativeOrder, derivativeOrder);

  int highestNonZeroCoefficientIndex = l;
  int numberNonZeroCoefficients = p_ - derivativeOrder + 1;
  int firstNonZeroIndex = highestNonZeroCoefficientIndex - numberNonZeroCoefficients + 1;
  const Eigen::VectorXd nonZeroCoefficients = allCoefficients.segment(firstNonZeroIndex, numberNonZeroCoefficients);

  Coefficients c{controlPointCount(), firstNonZeroIndex, nonZeroCoefficients};
  return c;
}

Eigen::VectorXd BSpline::front() const {
  return getControlPointMatrix().row(0).transpose();
}

Eigen::VectorXd BSpline::back() const {
  return getControlPointMatrix().row(controlPointCount() - 1).transpose();
}

Eigen::VectorXd BSpline::evaluate(double u, int derivativeOrder) const {
  if (derivativeOrder > p_) {
    return Eigen::VectorXd::Zero(dim_);
  }

  auto coefficients = calculateBSplineCoefficients(u, derivativeOrder);

  const auto& controlPoints = getControlPointMatrix(derivativeOrder);
  const auto& relevantControlPoints =
      controlPoints.middleRows(coefficients.firstNonZeroIndex(), coefficients.nonZeroCount());

  return relevantControlPoints.transpose() * coefficients.nonZeroCoefficients();
}

Eigen::VectorXd BSpline::deBoorCoefficient(double u, int i, int p, int k) const {
  if (p == 0) {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(controlPointCount());
    result(i) = 1;
    return result;
  }

  double alpha = (u - Uk_[k](i)) / (Uk_[k](i + (p_ - k) + 1 - p) - Uk_[k](i));
  Eigen::VectorXd result = (1 - alpha) * deBoorCoefficient(u, i - 1, p - 1, k) + alpha * deBoorCoefficient(u, i, p - 1, k);
  return result;
}

int BSpline::findIdxOfLowerOrEqualDomainKnot(double u, int k) const {
  assert(k <= p_);

  // start at the beginning of the domain [u_{p-k}^{(k)},u_{n+1-k}^{(k)}]
  int i = p_ - k;
  while ((u >= Uk_[k](i + 1)) && ((i + 1) < (n_ + 1 - k))) {
    ++i;
  }
  return i;
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
