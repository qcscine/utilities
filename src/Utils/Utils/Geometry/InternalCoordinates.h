/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_INTERNALCOORDINATES_H_
#define UTILS_INTERNALCOORDINATES_H_

#include "Utils/Geometry/AtomCollection.h"
#include <Eigen/Core>
#include <memory>

namespace Scine {
namespace Utils {

/**
 * @brief This class provides the option to convert from Cartesian into internal
 *        coordinates.
 *
 * As default internal redundant coordinates will be used, as a backup a set of
 * internal coordinates generated by only removing translation and rotational
 * degrees of freedom can be chosen.
 *
 * The redundant internal coordinates are generated using 'libirc':
 * https://github.com/RMeli/irc/
 */
class InternalCoordinates {
 public:
  /**
   * @brief Constructor
   * @param atoms The atoms for which the initial internal coordinates shall be
   *              generated.
   * @param rotTransOnly If true will only remove rotation and translation,
   *                     without actually generating internal coordinates from
   *                     bonds, angles and dihedrals.
   */
  InternalCoordinates(const AtomCollection& atoms, bool rotTransOnly = false);
  /// @brief Destructtor.
  ~InternalCoordinates();
  /**
   * @brief Back transforms internal coordinates into Cartesian ones.
   *
   * This procedure is iterative.
   *
   * @param internals The internal representation of the coordinates.
   * @param maxIters The maximum number of iterations to be tried for the back
   *                 transformation.
   * @param tolerance The tolerance (RMS of the internals) at which the back
   *                  transformation is considered converged.
   * @return PositionCollection The Cartesian representation of the coordinates.
   */
  PositionCollection coordinatesToCartesian(const Eigen::VectorXd& internals, unsigned int maxIters = 25,
                                            double tolerance = 1e-6) const;
  /**
   * @brief Transforms Cartesian coordinates into internal coordinates.
   * @param cartesian The Cartesian representation of the coordinates.
   * @return Eigen::VectorXd The internal representation of the coordinates.
   */
  Eigen::VectorXd coordinatesToInternal(const PositionCollection& cartesian) const;
  /**
   * @brief Transforms Cartesian gradients into internal gradients.
   *
   * This step does include the projection of the gradient into a valid one, for
   * the case of internal redundant coordinates.
   *
   * @param cartesian The Cartesian representation of the gradients.
   * @return Eigen::VectorXd The internal representation of the gradients.
   */
  Eigen::VectorXd gradientsToInternal(const GradientCollection& cartesian) const;
  /**
   * @brief Project the Hessian inverse such that it represents valid changes in
   *        internal coordinates.
   *
   * If only translation and rotation were removed this call does nothing to
   * the given matrix. It is thus safe to call this function in any case.
   *
   * Note: this projection is not a projection into internal coordinates from
   * Cartesian coordinates.
   *
   * @param inverse The Hessian inverse.
   * @return Eigen::MatrixXd The projected Hessian inverse.
   */
  Eigen::MatrixXd projectHessianInverse(const Eigen::MatrixXd& inverse) const;
  /**
   * @brief Project the Hessian such that it represents valid changes in
   *        internal coordinates.
   *
   * If only translation and rotation were removed this call does nothing to
   * the given matrix. It is thus safe to call this function in any case.
   *
   * Note: this projection is not a projection into internal coordinates from
   * Cartesian coordinates.
   *
   * @param hessian The Hessian.
   * @return Eigen::MatrixXd The projected Hessian.
   */
  Eigen::MatrixXd projectHessian(const Eigen::MatrixXd& hessian) const;
  /**
   * @brief Generates a guess for an inverse Hessian, the guess is a diagonal
   *        Matrix.
   *
   * If the internal representation corresponds only to the removal of
   * translation and rotation modes, the guess will be the identity.
   *
   * @return Eigen::MatrixXd The initial guess for the inverse Hessian.
   */
  Eigen::MatrixXd inverseHessianGuess() const;
  /**
   * @brief Generates a guess for a Hessian, the guess is a diagonal Matrix.
   *
   * If the internal representation corresponds only to the removal of
   * translation and rotation modes, the guess will be the identity.
   *
   * @return Eigen::MatrixXd The initial guess for the Hessian.
   */
  Eigen::MatrixXd hessianGuess() const;

 private:
  // Forward declaration for pointer to implementation idiom
  struct Impl;
  // Pointer to implementation in order to keep the YAML-cpp dependency private
  std::unique_ptr<Impl> _pImpl;
  // The initial cartesian coordinates needed for the back transformation
  mutable Eigen::VectorXd _oldCartesian;
  mutable Eigen::VectorXd _oldInternal;
};

/**
 * @brief Exception thrown when internal coordinates break down
 *
 * This exception is thrown if the interconversion between internal
 * and cartesian coordinates breaks down.
 */
class InternalCoordinatesException : public std::runtime_error {
 public:
  explicit InternalCoordinatesException()
    : std::runtime_error("Internal coordinates broke down, please try cartesians.") {
  }
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_INTERNALCOORDINATES_H_
