/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_GEOMETRYUTILITIES_H_
#define UTILS_GEOMETRYUTILITIES_H_

#include "Utils/Typenames.h"
#include <Eigen/Core>

namespace Scine {
namespace Utils {

class AtomCollection;
class Atom;

/// @brief Functionalities working with an entire geometry (PositionCollection).
namespace Geometry {

/**
 * @brief Translates a set of postions by a given displacement.
 * @param positions The original positions.
 * @param translation The displacement to be added.
 * @return PositionCollection Returns the translated positions.
 */
PositionCollection translatePositions(const PositionCollection& positions, const Eigen::Ref<Eigen::RowVector3d>& translation);
/**
 * @brief Translates a set of positions by a given displacement.
 *
 * The translation happens in-place.
 *
 * @param positions The initial positions, will be transformed in-place.
 * @param translation The displacement to be added.
 */
void translatePositionsInPlace(PositionCollection& positions, const Eigen::Ref<Eigen::RowVector3d>& translation);
/**
 * @brief Get the index of closest position (atom) to a given position in space.
 * @param positions A set of positions to be traversed.
 * @param targetPosition The target position.
 * @param distanceConsideredZero Squared distance between two positions resulting in them being considered equal
 *                               and therefore skipping this position as it is already present in the given
 *                               PositionCollection. The default is set to a negative value, so that all positions
 *                               are considered as possible candidates for the closest one.
 * @return int The index of the closest atom to the given target.
 */
int getIndexOfClosestAtom(const PositionCollection& positions, const Position& targetPosition,
                          double squaredDistanceConsideredZero = -1.0);
/**
 * @brief Get the index of a given atom in a given molecular structure.
 * @param structure The structure to be traversed.
 * @param atom The target atom that is tried to be found in the given structure.
 * @param distanceConsideredZero Squared distance between two atoms resulting in them being considered equal.
 * @return int The index of the given atom in the given molecular structure.
 * @throws std::runtime_error Function throws error if the atom is not found in the given structure.
 */
int getIndexOfAtomInStructure(const AtomCollection& structure, const Atom& atom, double squaredDistanceConsideredZero = 1e-4);
/**
 * @brief Transforms a 3N-dimensional vector {x0, y0, z0, x1, y1, z1, ...} to a Nx3 matrix.
 * @param v The positions in vector form.
 * @return Eigen::MatrixXd Returns the final Nx3 matrix.
 */
Eigen::MatrixXd positionVectorToMatrix(const Eigen::VectorXd& v);
/**
 * @brief Transforms a Nx3 matrix to a 3N-dimensional vector {x0, y0, z0, x1, y1, z1, ...}.
 * @param m The positions in matrix form.
 * @return Eigen::VectorXd Returns the final vector.
 */
Eigen::VectorXd positionMatrixToVector(const Eigen::MatrixXd& m);
/**
 * @brief Rotate and translate positions so that it is as close as possible to referencePositions.
 * @param reference The reference positions.
 * @param positions The positions to be aligned, will be transformed in place.
 */
void alignPositions(const PositionCollection& reference, PositionCollection& positions);
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
/**
 * @brief Calculated the cartesian modes corresponding to translations and rotations of the entire system.
 *
 * @param positions The positions of all atoms.
 * @param elements  The ElemenetTypes of all atoms.
 * @return Eigen::MatrixXd The rotation and translation modes.
 *                         Translation modes (x,y,z) first-third column,
 *                         rotation modes 3-final column (depending on the geometry)
 */
Eigen::MatrixXd calculateTranslationAndRotationModes(const PositionCollection& positions, const ElementTypeCollection& elements);

/**
 * @brief Generates the matrix removing rotation and translation modes from the given geometry if applied.
 *
 * @param positions The positions of all atoms.
 * @param elements  The ElemenetTypes of all atoms.
 * @param massWeighted True if Hessian to be transformed is mass-weighted
 * @return Eigen::MatrixXd The transformation matrix (applied as X^T*H*X to the Hessian).
 */
Eigen::MatrixXd calculateRotTransFreeTransformMatrix(const PositionCollection& positions,
                                                     const ElementTypeCollection& elements, bool massWeighted = false);

} /* namespace Geometry */
} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_GEOMETRYUTILITIES_H_
