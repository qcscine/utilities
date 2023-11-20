/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_BSPLINETOOLS_H
#define BSPLINES_BSPLINETOOLS_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {

namespace BSplines {

/*!
 * Contains methods that can be used by several B-spline classes.
 */
class BSplineTools {
 public:
  BSplineTools();

  static int findIdxOfLeftOrEqualDomainKnot(double u, int p, const Eigen::VectorXd& U);
  static int findIdxOfLeftDomainKnot(double u, int p, const Eigen::VectorXd& U);
  static int findIdxOfRightDomainKnot(double u, int p, const Eigen::VectorXd& U);
  static int findIdxOfRightOrEqualDomainKnot(double u, int p, const Eigen::VectorXd& U);

  static double knotAverage(int i, int p, const Eigen::VectorXd& U);

  static void normalizeKnotVector(Eigen::VectorXd& knotVector);

  static Eigen::VectorXd normalizedKnotVector(const Eigen::VectorXd& knotVector);

  static double rescaledKnot(double knot, std::pair<double, double> oldLim, std::pair<double, double> newLim);

  static void rescaleKnotVector(Eigen::VectorXd& knotVector, std::pair<double, double> oldLim, std::pair<double, double> newLim);

  static Eigen::VectorXd rescaledKnotVector(const Eigen::VectorXd& knotVector, std::pair<double, double> oldLim,
                                            std::pair<double, double> newLim);

  /* finite difference matrix penalizing BSplines */
  static int differenceOperator(int i, int j, int kappa);
};

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_BSPLINETOOLS_H
