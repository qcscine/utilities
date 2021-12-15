/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_KNOTINSERTER_H
#define BSPLINES_KNOTINSERTER_H

namespace Scine {
namespace Utils {

namespace BSplines {
class BSpline;

/*!
 * Implementation of the knot insertion algorithm by Böhm.
 * doi: 10.1016/0010-4485(80)90154-2
 */
struct KnotInserter {
  static void insertKnotByReference(double uInsert, BSpline& bs);
  static BSpline insertKnotByCopy(double u, const BSpline& bs);
};

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_KNOTINSERTER_H
