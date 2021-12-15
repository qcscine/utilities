/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_GEOMETRY_TRANSFORMATIONS_H_
#define UTILS_GEOMETRY_TRANSFORMATIONS_H_

#include "Utils/MolecularTrajectory.h"
#include "Utils/Typenames.h"
#include <Eigen/Core>
#include <Eigen/Dense>

namespace Scine {
namespace Utils {

/// @brief Functionalities working with an entire geometry (PositionCollection).
namespace Geometry {

/// @brief Functionalities working with changing datatypes or eliminating rotations and translations of an entire
/// geometry.
namespace Transformations {

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
 * @param elements  The ElementTypes of all atoms.
 * @param massWeighted True if Hessian to be transformed is mass-weighted
 * @return Eigen::MatrixXd The transformation matrix (applied as X^T*H*X to the Hessian).
 */
Eigen::MatrixXd calculateRotTransFreeTransformMatrix(const PositionCollection& positions,
                                                     const ElementTypeCollection& elements, bool massWeighted = false);

/**
 * @brief Generates the matrix removing rotation, translation modes and the component along a gradient from the given
 * geometry if applied.
 *
 * @param positions The positions of all atoms.
 * @param elements  The ElemenetTypes of all atoms.
 * @param gradients The cartesian gradient against which to orthogonalize.
 * @param massWeighted True if Hessian to be transformed is mass-weighted
 * @return Eigen::MatrixXd The transformation matrix (applied as X^T*H*X to the Hessian).
 */
Eigen::MatrixXd calculateRotTransFreeTransformMatrix(const PositionCollection& positions, const ElementTypeCollection& elements,
                                                     const GradientCollection& gradients, bool massWeighted = false);

} /* namespace Transformations */
} /* namespace Geometry */
} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_GEOMETRY_TRANSFORMATIONS_H_
