/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_STRUCTURALCOMPLETION_H_
#define UTILS_STRUCTURALCOMPLETION_H_

#include "Utils/Constants.h"
#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @class StructuralCompletion StructuralCompletion.h
 * @brief Contains function to generate positions from other positions.
 *
 * F.i. to be applied to add hydrogen atoms to a carbon chain.
 */
class StructuralCompletion {
 public:
  /// @brief Static functions only.
  StructuralCompletion() = delete;
  /**
   * @brief Generates three missing position in a tetrahedron.
   *
   * Given a normalized vector v1 showing from the center of a tetrahedron to one of its corners, calculates three
   * other normalized vectors to three other possible corners. Application example: -C -> -CH3.
   *
   * @param v1 A vector connecting the center and a known corner.
   * @param v2 A vector to be populated, then connecting the center and a corner.
   * @param v3 A vector to be populated, then connecting the center and a corner.
   * @param v4 A vector to be populated, then connecting the center and a corner.
   */
  static void generate3TetrahedronCornersFrom1Other(const Eigen::Vector3d& v1, Eigen::Ref<Eigen::Vector3d> v2,
                                                    Eigen::Ref<Eigen::Vector3d> v3, Eigen::Ref<Eigen::Vector3d> v4);
  /**
   * @brief Generates two missing position in a tetrahedron.
   *
   * Given two normalized vectors v1 and v2, calculates two other normalized vectors in order to generate a tetrahedron
   * (as well as possible). Does not work if v1 and v2 are collinear. Application example: -C- -> -CH2-.
   *
   * @param v1 A vector connecting the center and a known corner.
   * @param v2 A vector connecting the center and a known corner.
   * @param v3 A vector to be populated, then connecting the center and a corner.
   * @param v4 A vector to be populated, then connecting the center and a corner.
   */
  static void generate2TetrahedronCornersFrom2Others(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
                                                     Eigen::Ref<Eigen::Vector3d> v3, Eigen::Ref<Eigen::Vector3d> v4);
  /**
   * @brief Generates one missing position in a tetrahedron.
   *
   * Given three normalized vectors, calculates another normalized vector in order to generate a tetrahedron (as well
   * as possible). Does not work if v1, v2 and v3 are in the same plane. Application example: -C-R2 -> -CH-R2.
   *
   * @param v1 A vector connecting the center and a known corner.
   * @param v2 A vector connecting the center and a known corner.
   * @param v3 A vector connecting the center and a known corner.
   * @param v4 A vector to be populated, then connecting the center and a corner.
   */
  static void generate1TetrahedronCornerFrom3Others(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
                                                    const Eigen::Vector3d& v3, Eigen::Ref<Eigen::Vector3d> v4);
  /**
   * @brief Generates two missing position in a triangle.
   *
   * Given a normalized vector v1 showing from the center of a triangle to one of its corners, calculates two other
   * normalized vectors to two other possible corners. Application example: -C -> -CH2.
   *
   * @param v1 A vector connecting the center and a known corner.
   * @param v2 A vector to be populated, then connecting the center and a corner.
   * @param v3 A vector to be populated, then connecting the center and a corner.
   */
  static void generate2TriangleCornersFrom1Other(const Eigen::Vector3d& v1, Eigen::Ref<Eigen::Vector3d> v2,
                                                 Eigen::Ref<Eigen::Vector3d> v3);
  /**
   * @brief Generates a missing position in a triangle.
   *
   * Given two normalized vectors, calculates another normalized vector in order to generate a triangle (as well as
   * possible). Does not work if v1 and v2 are collinear. Application example: -C- -> -CH-.
   *
   * @param v1 A vector connecting the center and a known corner.
   * @param v2 A vector connecting the center and a known corner.
   * @param v3 A vector to be populated, then connecting the center and a corner.
   */
  static void generate1TriangleCornerFrom2Others(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
                                                 Eigen::Ref<Eigen::Vector3d> v3);

  /// @brief The standard tetrahedral angle in rad.
  static constexpr double tetrahedronAngle = 109.4712 * Constants::rad_per_degree;
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_STRUCTURALCOMPLETION_H_
