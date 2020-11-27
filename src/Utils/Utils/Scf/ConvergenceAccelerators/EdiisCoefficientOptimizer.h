/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_EDIISCOEFFICIENTOPTIMIZER_H
#define UTILS_EDIISCOEFFICIENTOPTIMIZER_H

#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {

/*!
 * Calculates optimal coefficients to use in EDIIS.
 * It solves the optimization problem given the B matrix and energies
 * and calculates the coefficients with the Reduced Gradient Method.
 */

class EdiisCoefficientOptimizer {
 public:
  EdiisCoefficientOptimizer(const std::vector<double>& energies, Eigen::MatrixXd B);

  /*! Calculate the coefficients. */
  Eigen::VectorXd getCoefficients();

 private:
  const Eigen::MatrixXd B_;
  Eigen::VectorXd energies_;
  unsigned dimension_;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_EDIISCOEFFICIENTOPTIMIZER_H