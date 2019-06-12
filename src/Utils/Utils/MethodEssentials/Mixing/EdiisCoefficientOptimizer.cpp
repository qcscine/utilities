/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "EdiisCoefficientOptimizer.h"
#include "ECQPP.h"

namespace Scine {
namespace Utils {

EdiisCoefficientOptimizer::EdiisCoefficientOptimizer(const std::vector<double>& energies, Eigen::MatrixXd B)
  : B_(std::move(B)) {
  dimension_ = static_cast<unsigned>(B_.rows());

  energies_.resize(dimension_);
  for (unsigned i = 0; i < dimension_; ++i)
    energies_[i] = energies[i];
}

Eigen::VectorXd EdiisCoefficientOptimizer::getCoefficients() {
  /*
   * Used to employ ReducedGradientMethod::Solver here, but it was sometimes incredibly slow to converge (millions
   * iterations) and sometimes found non-optimal minima.
   */
  ECQPP solver(B_, energies_);
  auto solution = solver.calculateOptimalCoefficients();
  return solution;
}
} // namespace Utils
} // namespace Scine
