/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Splitter.h"
#include "BSpline.h"
#include "BSplineTools.h"

namespace Scine {
namespace Utils {

namespace BSplines {

std::pair<BSpline, BSpline> Splitter::split(const double u, const BSpline& bs, const std::pair<bool, bool> normalizeKnotVectors) {
  assert((u > 0.0) && (u < 1.0));

  int p = bs.getDegree();
  BSpline bsInserted;
  bsInserted = bs;

  int l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(u, p, bs.getKnotVector());
  int it = l;
  int multiplicity = 0;

  while (bs.getKnotVector()[it] == u) {
    multiplicity++;
    it++;
  }

  assert(multiplicity <= p && "Spline cannot be disconnected");

  for (int i = 0; i < p - multiplicity; ++i) {
    KnotInserter::insertKnotByReference(u, bsInserted);
  }

  if (multiplicity > 0) {
    Eigen::VectorXd U(bsInserted.getKnotVector());
    Eigen::MatrixXd P(bsInserted.getControlPointMatrix());

    // left b-spline curve
    Eigen::VectorXd Uleft(l + p + 1);
    Uleft.head(l + p) = U.head(l + p);
    Uleft(l + p) = u;

    if (normalizeKnotVectors.first) {
      BSplineTools::normalizeKnotVector(Uleft);
    }

    Eigen::MatrixXd Pleft(P.topRows(Uleft.size() - (p + 2) + 1));

    // right b-spline curve
    Eigen::VectorXd Uright(U.size() - l + 1); // length of U - the elements that occur before the split idx
    Uright(0) = u;
    Uright.tail(U.size() - (l)) = U.tail(U.size() - (l));

    if (normalizeKnotVectors.second) {
      BSplineTools::normalizeKnotVector(Uright);
    }

    Eigen::MatrixXd Pright(P.bottomRows(Uright.size() - (p + 2) + 1));

    return std::make_pair(BSpline(Uleft, Pleft, p), BSpline(Uright, Pright, p));
  }

  Eigen::VectorXd U(bsInserted.getKnotVector());
  Eigen::MatrixXd P(bsInserted.getControlPointMatrix());

  // left b-spline curve
  Eigen::VectorXd Uleft(l + p + 1 + 1);
  Uleft.head(l + p + 1) = U.head(l + p + 1);
  Uleft(l + p + 1) = u;

  if (normalizeKnotVectors.first) {
    BSplineTools::normalizeKnotVector(Uleft);
  }

  // Eigen::MatrixXd Pleft(P.topRows(l+1));
  Eigen::MatrixXd Pleft(P.topRows(Uleft.size() - (p + 2) + 1));

  // right b-spline curve
  Eigen::VectorXd Uright(U.size() - (l + 1) + 1); // length of U - the elements that occur before the split idx
  Uright(0) = u;
  Uright.tail(U.size() - (l + 1)) = U.tail(U.size() - (l + 1));

  if (normalizeKnotVectors.second) {
    BSplineTools::normalizeKnotVector(Uright);
  }

  // Eigen::MatrixXd Pright(P.bottomRows(P.rows()-l));
  Eigen::MatrixXd Pright(P.bottomRows(Uright.size() - (p + 2) + 1));

  return std::make_pair(BSpline(Uleft, Pleft, p), BSpline(Uright, Pright, p));
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
