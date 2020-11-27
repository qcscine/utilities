/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_EXCEPTIONS_H
#define BSPLINES_EXCEPTIONS_H

#include <stdexcept>
#include <string>

namespace Scine {
namespace Utils {

namespace BSplines {

/*!
 * Base class for exceptions related to the BSplines library.
 */
class Exception : public std::runtime_error {
 public:
  explicit Exception(const std::string& s) : std::runtime_error(s) {
  }
};

class UnavailableDerivativeException : public Exception {
 public:
  UnavailableDerivativeException(int requestedDerivativeOrder, int highestAvailableDerivative)
    : Exception("Accessing BSpline derivatives of order " + std::to_string(requestedDerivativeOrder) +
                ", while only derivatives up to order " + std::to_string(highestAvailableDerivative) + " have been calculated.") {
  }
};

class InvalidSplitDimension : public Exception {
 public:
  InvalidSplitDimension() : Exception("The dimension index given in DimensionManipulator::split is invalid") {
  }
};

class IncompatibleKnotVectorsForDimensionalMerge : public Exception {
 public:
  IncompatibleKnotVectorsForDimensionalMerge()
    : Exception("The two BSplines being merged along dimensions have incompatible knot vectors.") {
  }
};

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_EXCEPTIONS_H