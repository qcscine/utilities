/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_OVERLAPCALCULATOR_H
#define UTILS_OVERLAPCALCULATOR_H

#include <Eigen/Core>

namespace Scine {

namespace Utils {
enum class DerivativeOrder;
}

namespace Utils {
class MatrixWithDerivatives;

/*!
 * Abstract class for the calculation of the overlap matrix.
 */
class OverlapCalculator {
 public:
  virtual ~OverlapCalculator() = default;

  /*! Calculate the overlap matrix upd to highestRequiredOrder. */
  virtual void calculateOverlap(Utils::DerivativeOrder highestRequiredOrder) = 0;

  virtual const MatrixWithDerivatives& getOverlap() const = 0;

  virtual void reset() = 0;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_OVERLAPCALCULATOR_H
