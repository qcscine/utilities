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
#include <vector>

namespace Scine {
namespace Utils {

//! @brief enum class defining the restrictedness of the reference calculation.
enum class Reference { Unrestricted, Restricted };

/**
 * @struct Container for any property depending on restricted/unrestricted.
 */
template<Reference restrictedness, typename Type>
struct SpinAdaptedContainer {};

/**
 * @struct Specialization for restricted containers.
 * This struct guarantees access to a restricted member.
 */
template<typename Type>
struct SpinAdaptedContainer<Reference::Restricted, Type> {
  SpinAdaptedContainer() = default;
  explicit SpinAdaptedContainer(Type instance) {
    restricted = std::move(instance);
  }
  Type restricted;
};

/**
 * @struct Specialization for unrestricted containers.
 * This struct guarantees access to alpha and beta members.
 */
template<typename Type>
struct SpinAdaptedContainer<Reference::Unrestricted, Type> {
  SpinAdaptedContainer() = default;
  explicit SpinAdaptedContainer(Type instance) {
    alpha = instance;
    beta = std::move(instance);
  }
  SpinAdaptedContainer(Type a, Type b) {
    alpha = std::move(a);
    beta = std::move(b);
  }
  Type alpha;
  Type beta;
};
//! @brief Alias for a pair of eigenvalues and eigenvectors.
struct EigenContainer {
  Eigen::VectorXd eigenValues;
  Eigen::MatrixXd eigenVectors;
};
/**
 * @brief Data structure to store the results of an excited states calculation.
 */
struct ElectronicTransitionResult {
  EigenContainer eigenStates;
  Eigen::Matrix3Xd transitionDipoles;
  Eigen::VectorXd spinContamination;
};

/**
 * @brief Data structure to store the results of an excited states calculation with closed-shell reference.
 */
struct SpinAdaptedElectronicTransitionResult {
  std::shared_ptr<ElectronicTransitionResult> unrestricted;
  std::shared_ptr<ElectronicTransitionResult> singlet;
  std::shared_ptr<ElectronicTransitionResult> triplet;
  std::vector<std::string> transitionLabels;
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
