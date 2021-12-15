/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_SPLITTER_H
#define BSPLINES_SPLITTER_H

#include "KnotInserter.h"
#include <Eigen/Core>
#include <utility>

namespace Scine {
namespace Utils {

namespace BSplines {
class BSpline;

/*! Splits a B-spline curve at a give parameter.
 * */
struct Splitter {
  static std::pair<BSpline, BSpline> split(double u, const BSpline& bs, std::pair<bool, bool> normalizeKnotVectors);
};

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_SPLITTER_H
