/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_BSPLINE_H
#define BSPLINES_BSPLINE_H

#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {

namespace BSplines {
class Coefficients;

/*!
 * Contains the standard implementation of a basis-spline curve.
 * Evaluation and getter methods are defined so that knot vectors and control points of derivatives are calculated
 * when they are needed.
 */
class BSpline {
 public:
  BSpline();
  BSpline(Eigen::VectorXd knotVector, Eigen::MatrixXd controlPoints, int degree);

  Eigen::VectorXd calculateBSplineCoefficientVector(double u, int derivativeOrder = 0) const;
  Coefficients calculateBSplineCoefficients(double u, int derivativeOrder = 0) const;

  /*! Return the value at u = 0.0. */
  Eigen::VectorXd front() const;
  /*! Return the value at u = 1.0. */
  Eigen::VectorXd back() const;

  Eigen::VectorXd evaluate(double u, int derivativeOrder = 0) const;
  Eigen::VectorXd operator()(double u) const;

  Eigen::VectorXd evaluateNaive(double u, int derivativeOrder = 0) const;

  void reverse();
  BSpline reversed() const;

  int getDegree() const {
    return p_;
  }

  int getOrder() const {
    return p_ + 1;
  }

  int getDim() const {
    return dim_;
  }

  int controlPointCount() const {
    return static_cast<int>(getControlPointMatrix().rows());
  }

  int getHighestCalculatedDerivative() const {
    return highestCalculatedDerivative_;
  }

  BSpline getDerivativeBSpline(int derivativeOrder);
  const Eigen::VectorXd& getKnotVector() const;
  const Eigen::VectorXd& getKnotVectorDerivative(int derivativeOrder) const;
  const Eigen::MatrixXd& getControlPointMatrix() const;
  const Eigen::VectorXd& deriveAndGetKnotVector(int derivativeOrder);
  const Eigen::MatrixXd& deriveAndGetControlPointMatrix(int derivativeOrder);
  bool isClampedAndNormed() const;

 private:
  void calculateDerivatives(int highestDerivativeToCalculate) const; // TODO: private!
  const Eigen::MatrixXd& getControlPointMatrix(int k) const;
  Eigen::VectorXd deBoorCoefficient(double u, int i, int p, int k = 0) const;
  void deriveKnotVectors(int highestDerivativeToCalculate) const;
  Eigen::VectorXd deriveControlPoint(int i, int highestDerivativeToCalculate) const;
  void deriveControlPointMatrices(int highestDerivativeToCalculate) const;
  int findIdxOfLowerOrEqualDomainKnot(double u, int k) const;

  int p_{0}, n_{0};
  int dim_{0};

  mutable int highestCalculatedDerivative_ = 0;
  mutable std::vector<Eigen::VectorXd> Uk_; // Knot vector and its derivatives
  mutable std::vector<Eigen::MatrixXd> Pk_; // Control points and their derivatives
};

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_BSPLINE_H
