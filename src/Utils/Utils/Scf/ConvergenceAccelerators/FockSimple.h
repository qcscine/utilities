/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_FOCK_SIMPLE_H
#define UTILS_FOCK_SIMPLE_H

#include <Utils/Scf/MethodInterfaces/ScfModifier.h>
#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {

/**
 * @class Fock_Simple fock_simple.h
 *
 * @brief A simple, damping-based convergence accelerator.
 *
 * This class is used to accelerate SCF convergence by implementing a damping scheme, which builds a linear
 * combination of the Fock matrices of the current and the previous SCF iterations. The amount with which
 * the Fock matrix of the previous iteration is mixed into the current one can be set via the damping
 * parameter.
 */
class FockSimple : public ScfModifier {
 public:
  /**
   * @brief This function is called from the corresponding SCF method in each SCF iteration after the calculation of
   * the Fock matrix.
   */
  void onFockCalculated() override;

  /**
   * @brief: Function to initialize the Fock_Simple object.
   */
  // TODO: Check whether this function can be made private.
  void initialize() override;

  /**
   * @brief Function to set the damping value (default is 0.8)
   *
   * @param damping The damping value to be used.
   */
  void setDamping(double damping) {
    damping_ = damping;
  }

  /**
   * @brief: Function to add a given Fock matrix to the list of stored Fock matrices.
   *
   * @param F Fock matrix to add.
   */
  // TODO: Check if this function can be made private.
  void addMatrices(const Eigen::MatrixXd& F);

 private:
  // This flag makes sure that the Fock_Simple object is only initialized once.
  bool initialized_{false};
  // The damping value to be used.
  double damping_{0.8};
  // The vector fockMatrices_ stores the Fock matrices of the current and the previous SCF iterations.
  std::vector<Eigen::MatrixXd> fockMatrices_;
  // The index (either 0 or 1) tells which element of fockMatrices_ is the Fock matrix of the current SCF iteration.
  int index_;
  // The number of atomic orbitals.
  int nAOs_;

  /*
   * @brief This function creates a linear combination of the Fock matrices of the current and the previous SCF
   * iterations.
   *
   * @return Eigen::MatrixXd& The new Fock matrix to be used in the next SCF iteration.
   */
  const Eigen::MatrixXd& extrapolate();
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_FOCK_SIMPLE_H
