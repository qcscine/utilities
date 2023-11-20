/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry/StructuralCompletion.h"
#include <Eigen/Geometry>

namespace Scine {
namespace Utils {

void StructuralCompletion::generate3TetrahedronCornersFrom1Other(const Eigen::Vector3d& v1, Eigen::Ref<Eigen::Vector3d> v2,
                                                                 Eigen::Ref<Eigen::Vector3d> v3,
                                                                 Eigen::Ref<Eigen::Vector3d> v4) {
  /*
   * Generate a vector perpendicular to the given one and then perform rotation around this vector to generate v2, then
   * rotate v2 around v1 to generate v3 and v4.
   */

  // Get a perpendicular vector;
  Eigen::Vector3d perpendicularVector = v1.cross(Eigen::Vector3d(1, 0, 0));
  // Check that it was different enough from {1,0,0}; if not, create it from {0,1,0}
  if (perpendicularVector.squaredNorm() < 0.00001) {
    perpendicularVector = v1.cross(Eigen::Vector3d(0, 1, 0));
  }
  perpendicularVector.normalize();

  Eigen::AngleAxisd t(tetrahedronAngle, perpendicularVector);
  auto R = t.toRotationMatrix();
  v2 = R * v1;
  double ang2 = 120.0 * Constants::rad_per_degree;
  Eigen::AngleAxisd t2(ang2, v1);
  R = t2.toRotationMatrix();
  v3 = R * v2;
  v4 = R * v3;
}

void StructuralCompletion::generate2TetrahedronCornersFrom2Others(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
                                                                  Eigen::Ref<Eigen::Vector3d> v3,
                                                                  Eigen::Ref<Eigen::Vector3d> v4) {
  /*
   * Generate the vector a, which is in the middle between v1 and v2, and the vector b, which is perpendicular to a but
   * in the same plane as v1 and v2, and then rotate a around b to generate v3 and v4.
   */

  Eigen::Vector3d a = (v1 + v2).normalized();
  Eigen::Vector3d b = (v2 - v1).normalized();

  double angle = (360 * Constants::rad_per_degree - tetrahedronAngle) / 2;
  Eigen::AngleAxisd t(angle, b);
  auto R = t.toRotationMatrix();
  v3 = R * a;
  v4 = R * v3;
}

void StructuralCompletion::generate1TetrahedronCornerFrom3Others(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
                                                                 const Eigen::Vector3d& v3, Eigen::Ref<Eigen::Vector3d> v4) {
  Eigen::Vector3d res = -(v1 + v2 + v3);
  // check that the three vectors were non-planar enough
  if (res.squaredNorm() < 0.4 * 0.4) {
    res = v1.cross(v2);
  }

  v4 = res.normalized();
}

void StructuralCompletion::generate2TriangleCornersFrom1Other(const Eigen::Vector3d& v1, Eigen::Ref<Eigen::Vector3d> v2,
                                                              Eigen::Ref<Eigen::Vector3d> v3) {
  /*
   * Generate a vector perpendicular to the given one and then perform rotation around this vector to generate v2 and
   * v3.
   */

  // Get a perpendicular vector;
  auto perpendicularVector = v1.cross(Eigen::Vector3d(1, 0, 0));
  // Check that it was different enough from {1,0,0}; if not, create it from {0,1,0}
  if (perpendicularVector.squaredNorm() < 0.00001) {
    perpendicularVector = v1.cross(Eigen::Vector3d(0, 1, 0));
  }
  perpendicularVector.normalize();

  Eigen::AngleAxisd t(120 * Constants::rad_per_degree, perpendicularVector);
  auto R = t.toRotationMatrix();
  v2 = R * v1;
  v3 = R * v2;
}

void StructuralCompletion::generate1TriangleCornerFrom2Others(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
                                                              Eigen::Ref<Eigen::Vector3d> v3) {
  v3 = -(v1 + v2).normalized();
}

constexpr double StructuralCompletion::tetrahedronAngle;

} /* namespace Utils */
} /* namespace Scine */
