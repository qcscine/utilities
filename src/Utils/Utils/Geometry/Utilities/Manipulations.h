/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_GEOMETRY_MANIPULATIONS_H_
#define UTILS_GEOMETRY_MANIPULATIONS_H_

#include "Utils/Geometry/AtomCollection.h"
#include "Utils/MolecularTrajectory.h"
#include "Utils/Typenames.h"
#include <Eigen/Core>
#include <Eigen/Dense>

namespace Scine {
namespace Utils {

class MolecularTrajectory;

/// @brief Functionalities working with an entire geometry (PositionCollection).
namespace Geometry {

/// @brief Functionalities to manipulate an entire geometry.
namespace Manipulations {
/**
 * @brief Translates a set of positions by a given displacement.
 * @param positions The original positions.
 * @param translation The displacement to be added.
 * @return PositionCollection Returns the translated positions.
 */
PositionCollection translatePositions(const PositionCollection& positions, const Eigen::RowVector3d& translation);
/**
 * @brief Translates a set of positions by a given displacement.
 *
 * The translation happens in-place.
 *
 * @param positions The initial positions, will be transformed in-place.
 * @param translation The displacement to be added.
 */
void translatePositionsInPlace(PositionCollection& positions, const Eigen::RowVector3d& translation);
/**
 * @brief Rotate position collection such that startOrientation vector overlaps with endOrientation vector.
 * @param positions Position collection of interest.
 * @param startOrientation Vector of start orientation of positions.
 * @param endOrientation Vector of end orientation of positions.
 * @param pointOfRotation Point through which both orientations pass.
 * @return Position collection of rotated positions.
 */
PositionCollection rotatePositions(const PositionCollection& positions, const Eigen::RowVector3d& startOrientation,
                                   const Eigen::RowVector3d& endOrientation, const Eigen::RowVector3d& pointOfRotation);
/**
 * @brief Rotate position collection around given axis by given angle.
 * @param positions Position collection of interest.
 * @param rotAxisOrientation Vector of rotation axis. Will be normalized.
 * @param angle Rotation angle in radians.
 * @param pointOfRotation Origin of rotation axis.
 * @return Position collection of rotated positions.
 */
PositionCollection rotatePositions(const PositionCollection& positions, const Eigen::RowVector3d& rotAxisOrientation,
                                   double angle, const Eigen::RowVector3d& pointOfRotation);
/**
 * @brief Rotate position collection by given quaternion.
 *
 * The rotation happens in-place.
 *
 * @param pc Position collection of interest.
 * @param rotation Quaternion<double> containing the information about the rotation.
 * @param pointOfRotation Point around which pc will be rotated.
 */
void rotatePositionsInPlace(PositionCollection& pc, const Eigen::Quaterniond& rotation,
                            const Eigen::RowVector3d& pointOfRotation);
/**
 * @brief Randomly displaces all positions in all directions.
 *
 * @param positions The positions of interest.
 * @param maxDisplacement The maximum displacement for each coordinate.
 * @return PositionCollection of the randomly displaced positions
 */
PositionCollection randomDisplacement(const PositionCollection& positions, double maxDisplacement);
/**
 * @brief Randomly displaces all positions in all directions.
 *
 * The displacement happens in-place.
 *
 * @param positions The initial positions, which will be displaced in-place.
 * @param maxDisplacement The maximum displacement for each coordinate.
 */
void randomDisplacementInPlace(PositionCollection& positions, double maxDisplacement);
/**
 * @brief Randomly displaces positions in AtomCollection and saves them as MolecularTrajectory.
 *
 * @param atoms The atoms of interest.
 * @param numFrames The number of frames in the returned trajectory.
 * @param maxDisplacement The maximum displacement for each coordinate in each frame.
 * @return MolecularTrajectory of the randomly displaced atoms.
 */
MolecularTrajectory randomDisplacementTrajectory(const AtomCollection& atoms, unsigned int numFrames, double maxDisplacement);
/**
 * @brief Randomly displaces positions in AtomCollection and saves them as MolecularTrajectory.
 *
 * @param atoms The atoms of interest.
 * @param numFrames The number of frames in the returned trajectory.
 * @param maxDisplacement The maximum displacement for each coordinate in each frame.
 * @param seed Specify the seed for the pseudorandom number generator.
 * @return MolecularTrajectory of the randomly displaced atoms.
 */
MolecularTrajectory randomDisplacementTrajectory(const AtomCollection& atoms, unsigned int numFrames,
                                                 double maxDisplacement, unsigned int seed);
/**
 * @brief Rotate and translate positions so that it is as close as possible to referencePositions.
 * @param reference The reference positions.
 * @param positions The positions to be aligned, will be transformed in place.
 */
void alignPositions(const PositionCollection& reference, PositionCollection& positions);
/**
 * @brief Rotate and translate positions so that it is as close as possible to referencePositions.
 * @param reference The reference positions.
 * @param positions The positions to be aligned, will be transformed in place.
 * @param elements The elements composing the structure to allow a mass-weighting of the alignment.
 */
void alignPositions(const PositionCollection& reference, PositionCollection& positions, const ElementTypeCollection& elements);
/**
 * @brief Rotate and translate positions so that it is as close as possible to referencePositions with weights.
 * @param reference The reference positions.
 * @param positions The positions to be aligned, will be transformed in place.
 * @param weights The weigths with which each atom goes into the alignment.
 * @param rmsdVector The vector in which the rmsd of the alignment are written.
 */
void alignPositions(const PositionCollection& reference, PositionCollection& positions, const Eigen::VectorXd& weights,
                    Eigen::VectorXd& rmsdVector);

} /* namespace Manipulations */
} /* namespace Geometry */
} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_GEOMETRY_MANIPULATIONS_H_
