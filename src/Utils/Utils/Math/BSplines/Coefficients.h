/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_COEFFICIENTS_H
#define BSPLINES_COEFFICIENTS_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {

namespace BSplines {

/*!
 * Class for the B-Spline coefficients N for some given derivativeOrder and u value.
 * Since many of the coefficients are equal to zero, just the non-zero coefficients are saved.
 */
class Coefficients {
 public:
  explicit Coefficients(int coefficientCount, int firstNonZeroIndex, Eigen::VectorXd nonZeroValues);

  Eigen::VectorXd fullCoefficientVector() const;
  const Eigen::VectorXd& nonZeroCoefficients() const;

  int firstNonZeroIndex() const;
  int lastNonZeroIndex() const;
  int nonZeroCount() const;

  double get(int index) const;

 private:
  int coefficientCount_;
  int firstNonZero_;
  Eigen::VectorXd nonZeroCoefficients_;
};

inline Coefficients::Coefficients(int coefficientCount, int firstNonZeroIndex, Eigen::VectorXd nonZeroValues)
  : coefficientCount_(coefficientCount), firstNonZero_(firstNonZeroIndex), nonZeroCoefficients_(std::move(nonZeroValues)) {
}

inline Eigen::VectorXd Coefficients::fullCoefficientVector() const {
  Eigen::VectorXd coefficients = Eigen::VectorXd::Zero(coefficientCount_);
  coefficients.segment(firstNonZero_, nonZeroCount()) = nonZeroCoefficients_;
  return coefficients;
}

inline const Eigen::VectorXd& Coefficients::nonZeroCoefficients() const {
  return nonZeroCoefficients_;
}

inline int Coefficients::firstNonZeroIndex() const {
  return firstNonZero_;
}

inline int Coefficients::lastNonZeroIndex() const {
  return firstNonZero_ + nonZeroCount() - 1;
}

inline int Coefficients::nonZeroCount() const {
  return static_cast<int>(nonZeroCoefficients_.size());
}

inline double Coefficients::get(int index) const {
  if (index < firstNonZero_) {
    return 0.0;
  }
  if (index < firstNonZero_ + nonZeroCoefficients_.size()) {
    return nonZeroCoefficients_(index - firstNonZero_);
  }
  return 0;
}

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_COEFFICIENTS_H
