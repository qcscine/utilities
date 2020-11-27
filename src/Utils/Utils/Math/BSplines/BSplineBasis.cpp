/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "BSplineBasis.h"

namespace Scine {
namespace Utils {

namespace BSplines {

/*
 * p is the degree of this basis function of the kth degree spline - not of the 0th degree spline
 * */
double BSplineBasis::evaluate(int i, int pk, int numberOfSplineSegments, const Eigen::VectorXd& knotVector, double u) {
  assert((u >= knotVector(pk)) && (u <= knotVector(numberOfSplineSegments + 1)) &&
         "Parameter has to lie inside of the domain [u_{p},u_{n+1}]");

  double summand1 = 0, summand2 = 0;
  if (pk == 0) {
    /*! the conditional (i == n) &&  (u == Uk(n+1)) closes the last interval [u_{n},u_{n+1}] */
    if (((knotVector(i) <= u) && (u < knotVector(i + 1))) ||
        ((i == numberOfSplineSegments) && (u == knotVector(numberOfSplineSegments + 1))))
      return 1;
    else
      return 0;
  }
  else {
    if (knotVector(i + pk) == knotVector(i))
      summand1 = 0;
    else
      summand1 = (u - knotVector(i)) / (knotVector(i + pk) - knotVector(i)) *
                 evaluate(i, pk - 1, numberOfSplineSegments, knotVector, u);

    if (knotVector(i + pk + 1) == knotVector(i + 1))
      summand2 = 0;
    else
      summand2 = (knotVector(i + pk + 1) - u) / (knotVector(i + pk + 1) - knotVector(i + 1)) *
                 evaluate(i + 1, pk - 1, numberOfSplineSegments, knotVector, u);

    return summand1 + summand2;
  }
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
