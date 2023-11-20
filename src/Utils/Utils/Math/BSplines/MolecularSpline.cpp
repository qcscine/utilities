/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MolecularSpline.h"
#include <Utils/Geometry/AtomCollection.h>

namespace Scine {
namespace Utils {

namespace BSplines {

Utils::PositionCollection MolecularSpline::getPositions(double u) const {
  Eigen::VectorXd v = spline_.evaluate(u);
  Utils::PositionCollection positions = Eigen::Map<const Utils::PositionCollection>(v.data(), v.size() / 3, 3);
  return positions;
}

Utils::AtomCollection MolecularSpline::at(double u) const {
  return Utils::AtomCollection(getElements(), getPositions(u));
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
