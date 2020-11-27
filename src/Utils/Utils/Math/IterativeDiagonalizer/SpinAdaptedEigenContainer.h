/**
 * @file EigenPairs.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EIGENPAIRS_H
#define UTILS_EIGENPAIRS_H

#include <Eigen/Core>
#include <memory>

namespace Scine {
namespace Utils {

//! @brief Alias for a pair of eigenvalues and eigenvectors.
using EigenContainer = std::pair<Eigen::VectorXd, Eigen::MatrixXd>;
/**
 * @brief Data structure to store the results of an excited states calculation.
 */
struct ElectronicTransitionResult {
  EigenContainer eigenStates;
  Eigen::Matrix3Xd transitionDipoles;
};

/**
 * @brief Data structure to store the results of an excited states calculation with closed-shell reference.
 */
struct SpinAdaptedElectronicTransitionResult {
  std::shared_ptr<ElectronicTransitionResult> singlet;
  std::shared_ptr<ElectronicTransitionResult> triplet;
};

/**
 * @brief Coupling matrix needed by TD-DFTB and TD-DFT
 */
struct CouplingMatrix {
  std::shared_ptr<Eigen::MatrixXd> singletK;
  std::shared_ptr<Eigen::MatrixXd> tripletK;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_EIGENPAIRS_H
