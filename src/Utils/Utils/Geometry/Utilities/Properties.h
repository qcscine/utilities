/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_GEOMETRY_PROPERTIES_H_
#define UTILS_GEOMETRY_PROPERTIES_H_

#include "Utils/Geometry/AtomCollection.h"
#include <Eigen/Core>
#include <Eigen/Dense>

namespace Scine {
namespace Utils {

/// @brief Functionalities working with an entire geometry (PositionCollection).
namespace Geometry {

/// @brief Functionalities giving properties of an entire geometry.
namespace Properties {

/**
 * @class PrincipalMomentsOfInertia Geometry.h
 * @brief The principal moments of inertia stored.
 */
struct PrincipalMomentsOfInertia {
  Eigen::Vector3d eigenvalues;
  Eigen::Matrix3d eigenvectors;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
};
/**
 * @brief Get a vector of all masses (in a.u.).
 * @param elements A collection of elements.
 * @return std::vector<double> Returns the masses listed in a vector.
 */
std::vector<double> getMasses(const ElementTypeCollection& elements);
/**
 * @brief Get the center of mass.
 * @param positions The positions.
 * @param masses The masses (sorted according to the positions).
 * @return Position Returns the center of mass (COM).
 */
Position getCenterOfMass(const PositionCollection& positions, const std::vector<double>& masses);
/**
 * @brief Get the center of mass.
 * @param structure The structure (positions and masses are relevant).
 * @return Position Returns the center of mass (COM).
 */
Position getCenterOfMass(const AtomCollection& structure);
/**
 * @brief Get the average Position.
 *
 * (The average Position is identical to the center of mass if all masses are identical)
 *
 * @param positions The positions.
 * @return Position Returns the average position.
 */
Position getAveragePosition(const PositionCollection& positions);
/**
 * @brief Calculates the inertia tensor.
 * @param positions The positions.
 * @param masses The masses (sorted according to the positions).
 * @param centerOfMass The center of mass.
 * @return Eigen::Matrix3d Returns the inertia tensor.
 */
Eigen::Matrix3d calculateInertiaTensor(const PositionCollection& positions, const std::vector<double>& masses,
                                       const Position& centerOfMass);
/**
 * @brief Calculates the principal moments of inertia.
 * @param positions The positions.
 * @param masses The masses (sorted according to the positions).
 * @param centerOfMass The center of mass.
 * @return PrincipalMomentsOfInertia Returns the principal moments of inertia.
 */
PrincipalMomentsOfInertia calculatePrincipalMoments(const PositionCollection& positions,
                                                    const std::vector<double>& masses, const Position& centerOfMass);

bool operator==(const PositionCollection& lhs, const PositionCollection& rhs);
bool operator!=(const PositionCollection& lhs, const PositionCollection& rhs);

} /* namespace Properties */
} /* namespace Geometry */
} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_GEOMETRY_PROPERTIES_H_
