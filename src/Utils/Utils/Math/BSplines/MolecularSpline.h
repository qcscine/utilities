/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_MOLECULARSPLINE_H
#define BSPLINES_MOLECULARSPLINE_H

#include "BSpline.h"
#include <Utils/Typenames.h>

namespace Scine {

namespace Utils {
class AtomCollection;
} // namespace Utils

namespace Utils {
namespace BSplines {

/*!
 * Class for a molecular trajectory saved as a B-spline.
 * In principle, just a struct of a BSpline and a ElementTypeCollection.
 */
class MolecularSpline {
 public:
  MolecularSpline() = default;
  MolecularSpline(Utils::ElementTypeCollection elements, BSpline spline);

  const Utils::ElementTypeCollection& getElements() const;
  BSpline& getBSpline();
  const BSpline& getBSpline() const;

  /*! Get the positions at a given coordinate along the spline. */
  Utils::PositionCollection getPositions(double u) const;
  /*! Generate an AtomCollection at a given coordinate along the spline. */
  Utils::AtomCollection at(double u) const;

 private:
  Utils::ElementTypeCollection elements_;
  BSpline spline_{};
};

inline MolecularSpline::MolecularSpline(Utils::ElementTypeCollection elements, BSpline spline)
  : elements_(std::move(elements)), spline_(std::move(spline)) {
}

inline const Utils::ElementTypeCollection& MolecularSpline::getElements() const {
  return elements_;
}

inline BSpline& MolecularSpline::getBSpline() {
  return spline_;
}

inline const BSpline& MolecularSpline::getBSpline() const {
  return spline_;
}

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_MOLECULARSPLINE_H