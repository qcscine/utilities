/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_CONVERGENCECHECKER_H
#define UTILS_CONVERGENCECHECKER_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {

class ScfMethod;

/*!
 * Class for a density matrix based convergence checker.
 * It looks at the RMS change in density matrix elements.
 */

class ConvergenceChecker {
 public:
  explicit ConvergenceChecker(const ScfMethod& method);

  /*! Set the threshold under which the convergence criterion is satisfied.
      Corresponds to the rmsd in the density matrix change. */
  void setThreshold(double v) {
    threshold_ = v;
  }
  /*! Get the currently applied threshold. */
  double getThreshold() const {
    return threshold_;
  }

  /*! Gets the current density matrix from the assigned SCF method.
      Performs no check of convergence and is useful if the starting density matrix in the SCF cycle is not the last one
     of the previous cycle. */
  void updateDensityMatrix();

  /*! Check if the convergence criterion is met.
      Implicitly calls updateDensityMatrix(). */
  void checkConvergence();

  bool isConverged() const {
    return converged_;
  }

 private:
  double calculateRMSD();

  const ScfMethod& m_;
  bool converged_;
  double threshold_;
  Eigen::MatrixXd oldMatrix_, newMatrix_;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_CONVERGENCECHECKER_H