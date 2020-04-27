/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "KnotInserter.h"
#include "BSpline.h"
#include "BSplineTools.h"

namespace Scine {
namespace Utils {

namespace BSplines {

KnotInserter::KnotInserter() = default;

/*!
 * Boehm's algorithm for knot insertion
 * the knot must be >0 and <1 */
void KnotInserter::insertKnotByReference(const double uInsert, BSpline& bs) const {
  assert((uInsert > 0.0) && (uInsert < 1.0));

  int p = bs.getDegree();
  int dim = bs.getDim();

  Eigen::VectorXd Uold = bs.getKnotVector();
  Eigen::MatrixXd Pold = bs.getControlPointMatrix();
  auto numberOfControlPoints = static_cast<int>(Pold.rows());
  auto numberOfKnots = static_cast<int>(Uold.size());

  // find knot span
  int l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(uInsert, p, Uold);

  // compute new control points. Only the control points l-p+1 to l change
  Eigen::MatrixXd Pnew(numberOfControlPoints + 1, dim);
  Pnew.topRows(l - p + 1) = Pold.topRows(l - p + 1);
  Pnew.bottomRows(numberOfControlPoints - l) = Pold.bottomRows(numberOfControlPoints - l);

  double alpha;
  for (int i = l - p + 1; i <= l; ++i) {
    alpha = (uInsert - Uold(i)) / (Uold(i + p) - Uold(i));
    Pnew.row(i) = (1 - alpha) * Pold.row(i - 1) + alpha * Pold.row(i);
  }

  // insert u into knot vector
  Eigen::VectorXd Unew(numberOfKnots + 1);
  Unew.head(l + 1) = Uold.head(l + 1);
  Unew(l + 1) = uInsert;
  Unew.tail(numberOfKnots - l - 1) = Uold.tail(numberOfKnots - l - 1);

  bs = BSpline(Unew, Pnew, p);
}

BSpline KnotInserter::insertKnotByCopy(const double u, const BSpline& bs) const {
  auto bsCopy = bs;
  KnotInserter::insertKnotByReference(u, bsCopy);
  return bsCopy;
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
