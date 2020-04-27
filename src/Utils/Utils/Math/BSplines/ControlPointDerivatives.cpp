/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ControlPointDerivatives.h"
#include "BSpline.h"

namespace Scine {
namespace Utils {

namespace BSplines {

Eigen::VectorXd ControlPointDerivatives::oneDerivative(const BSpline& spline, double u, int controlPointIndex,
                                                       int curveDerivative) {
  // The derivative for controlPointIndex is actually as if the control point matrix is zero everywhere but for
  // the considered control point.
  auto& cp = spline.getControlPointMatrix();

  Eigen::MatrixXd newCpMatrix = Eigen::MatrixXd::Zero(cp.rows(), cp.cols());
  newCpMatrix.row(controlPointIndex).setOnes();

  auto cpBSpline = BSpline{spline.getKnotVector(), newCpMatrix, spline.getDegree()};

  return cpBSpline.evaluate(u, curveDerivative);
}

Eigen::MatrixXd ControlPointDerivatives::allDerivatives(const BSpline& spline, double u, int curveDerivative) {
  Eigen::MatrixXd cpDerivatives(spline.controlPointCount(), spline.getDim());
  for (int i = 0; i < spline.controlPointCount(); ++i) {
    Eigen::VectorXd der = oneDerivative(spline, u, i, curveDerivative);
    cpDerivatives.row(i) = der.transpose();
  }
  return cpDerivatives;
}

Eigen::VectorXd ControlPointDerivatives::curveDerivatives(const BSpline& spline, double u, int controlPointIndex) {
  return oneDerivative(spline, u, controlPointIndex, 0);
}

Eigen::MatrixXd ControlPointDerivatives::curveDerivatives(const BSpline& spline, double u) {
  return allDerivatives(spline, u, 0);
}

Eigen::VectorXd ControlPointDerivatives::firstOrderCurveDerivatives(const BSpline& spline, double u, int controlPointIndex) {
  return oneDerivative(spline, u, controlPointIndex, 1);
}

Eigen::MatrixXd ControlPointDerivatives::firstOrderCurveDerivatives(const BSpline& spline, double u) {
  return allDerivatives(spline, u, 1);
}

Eigen::VectorXd ControlPointDerivatives::secondOrderCurveDerivatives(const BSpline& spline, double u, int controlPointIndex) {
  return oneDerivative(spline, u, controlPointIndex, 2);
}

Eigen::MatrixXd ControlPointDerivatives::secondOrderCurveDerivatives(const BSpline& spline, double u) {
  return allDerivatives(spline, u, 2);
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
