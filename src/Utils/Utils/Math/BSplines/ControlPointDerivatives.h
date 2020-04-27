/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_CONTROLPOINTDERIVATIVES_H
#define BSPLINES_CONTROLPOINTDERIVATIVES_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {

namespace BSplines {
class BSpline;

/*!
 * Calculation of the derivatives of BSpline curves C (or of dC/du, d2C/du2, ...) with respect
 * to the control points.
 *
 * Not only can the derivatives of the curve itself with respect to the control points be calculated, but also
 * of dC/du, d2C/du2, etc. The degree of the curve to derive is given by "curveDerivative" in some functions below.
 *
 * The functions exist in two variants: one to calculate the derivatives of the curve with respect to one control point,
 * or with respect to all the control points.
 *
 * Note that the derivative of some control point coordinate i is non-zero only for the coordinate i of the curve.
 */
namespace ControlPointDerivatives {

Eigen::VectorXd oneDerivative(const BSpline& spline, double u, int controlPointIndex, int curveDerivative);
Eigen::MatrixXd allDerivatives(const BSpline& spline, double u, int curveDerivative);

/*! Get the derivatives of the curve at a given u with respect to some control point coordinates.
 * returns a vector of dCk / dPk for k in [0, nDimensions). NB: all other dCi / dPj are zero. */
Eigen::VectorXd curveDerivatives(const BSpline& spline, double u, int controlPointIndex);

/*! Get the derivatives of the curve at a given u with respect to all control point coordinates.
 * returns a matrix of dCk / dP^A_k for k in [0, nDimensions), and A being a given control point.
 * NB: all other dCi / dPj are zero. */
Eigen::MatrixXd curveDerivatives(const BSpline& spline, double u);

/*! Get the derivatives of the curve tangent at a given u with respect to some control point coordinates.
 * returns a vector of dC'k / dPk for k in [0, nDimensions). NB: all other dC'i / dPj are zero. */
Eigen::VectorXd firstOrderCurveDerivatives(const BSpline& spline, double u, int controlPointIndex);

/*! Get the derivatives of the curve tangent at a given u with respect to all control point coordinates.
 * returns a matrix of dC'k / dP^A_k for k in [0, nDimensions), and A being a given control point.
 * NB: all other dC'i / dPj are zero. */
Eigen::MatrixXd firstOrderCurveDerivatives(const BSpline& spline, double u);

/*! Get the derivatives of the curve tangent at a given u with respect to some control point coordinates.
 * returns a vector of dC'k / dPk for k in [0, nDimensions). NB: all other dC'i / dPj are zero. */
Eigen::VectorXd secondOrderCurveDerivatives(const BSpline& spline, double u, int controlPointIndex);

/*! Get the derivatives of the curve tangent at a given u with respect to all control point coordinates.
 * returns a matrix of dC'k / dP^A_k for k in [0, nDimensions), and A being a given control point.
 * NB: all other dC'i / dPj are zero. */
Eigen::MatrixXd secondOrderCurveDerivatives(const BSpline& spline, double u);

} // namespace ControlPointDerivatives

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_CONTROLPOINTDERIVATIVES_H