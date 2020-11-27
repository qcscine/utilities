/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Math/AutomaticDifferentiation/VectorDerivatives3D.h>
#include <gmock/gmock.h>
#include <Eigen/Geometry>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

using namespace AutomaticDifferentiation;

/**
 * @class AVectorDerivatives3DTest VectorDerivatives3DTest.cpp
 * @brief Comprises tests for the class Scine::Utils::VectorDerivatives3D.
 * @test
 */
class AVectorDerivatives3DTest : public Test {
 public:
  double x = 0.231;
  double y = -0.2922;
  double z = 0.787;
  Eigen::Vector3d v;

 protected:
  void SetUp() override {
    v = Eigen::Vector3d(x, y, z);
  }
};

TEST_F(AVectorDerivatives3DTest, IsCorrectForDotProductWithItself) {
  auto dv = VectorDerivatives3D::spatialVectorHessian3D(v);

  AutomaticDifferentiation::Second3D dotProduct = dv.dot(dv);

  ASSERT_THAT(dotProduct.value(), DoubleEq(x * x + y * y + z * z));
  ASSERT_THAT(dotProduct.dx(), DoubleEq(2 * x));
  ASSERT_THAT(dotProduct.dy(), DoubleEq(2 * y));
  ASSERT_THAT(dotProduct.dz(), DoubleEq(2 * z));
  ASSERT_THAT(dotProduct.XX(), DoubleEq(2));
  ASSERT_THAT(dotProduct.YY(), DoubleEq(2));
  ASSERT_THAT(dotProduct.ZZ(), DoubleEq(2));
  ASSERT_THAT(dotProduct.XY(), DoubleEq(0));
  ASSERT_THAT(dotProduct.XZ(), DoubleEq(0));
  ASSERT_THAT(dotProduct.YZ(), DoubleEq(0));
}

TEST_F(AVectorDerivatives3DTest, IsCorrectForDotProductWithVector3d) {
  auto dv = VectorDerivatives3D::spatialVectorHessian3D(v);
  double x2 = -3.2, y2 = 0.21, z2 = 5.121;
  Eigen::Vector3d v2(x2, y2, z2);

  AutomaticDifferentiation::Second3D dotProduct = dv.dot(v2);

  ASSERT_THAT(dotProduct.value(), DoubleEq(x * x2 + y * y2 + z * z2));
  ASSERT_THAT(dotProduct.dx(), DoubleEq(x2));
  ASSERT_THAT(dotProduct.dy(), DoubleEq(y2));
  ASSERT_THAT(dotProduct.dz(), DoubleEq(z2));
  ASSERT_THAT(dotProduct.XX(), DoubleEq(0));
  ASSERT_THAT(dotProduct.YY(), DoubleEq(0));
  ASSERT_THAT(dotProduct.ZZ(), DoubleEq(0));
  ASSERT_THAT(dotProduct.XY(), DoubleEq(0));
  ASSERT_THAT(dotProduct.XZ(), DoubleEq(0));
  ASSERT_THAT(dotProduct.YZ(), DoubleEq(0));
}

TEST_F(AVectorDerivatives3DTest, IsCorrectForCrossProduct) {
  double x2 = -3.2, y2 = 0.21, z2 = 5.121;
  Eigen::Vector3d v2(x2, y2, z2);
  double factor = 3.121;

  auto dv = VectorDerivatives3D::spatialVectorHessian3D(v);
  auto dv2 = VectorDerivatives3D::spatialVectorHessian3D(v2) * factor;

  auto crossProduct = dv.cross(dv2);
  Eigen::Vector3d expectedValue = v.cross(v2 * factor);

  // Check values (compare with Eigen::cross):
  ASSERT_THAT(crossProduct.x().value(), DoubleEq(expectedValue.x()));
  ASSERT_THAT(crossProduct.y().value(), DoubleEq(expectedValue.y()));
  ASSERT_THAT(crossProduct.z().value(), DoubleEq(expectedValue.z()));

  // Check only x component derivatives
  auto xDer = crossProduct.x();
  ASSERT_THAT(xDer.dx(), DoubleEq(0));
  ASSERT_THAT(xDer.dy(), DoubleEq(factor * z2 - factor * z));
  ASSERT_THAT(xDer.dz(), DoubleEq(factor * y - factor * y2));
  ASSERT_THAT(xDer.XX(), DoubleEq(0));
  ASSERT_THAT(xDer.YY(), DoubleEq(0));
  ASSERT_THAT(xDer.ZZ(), DoubleEq(0));
  ASSERT_THAT(xDer.XY(), DoubleEq(0));
  ASSERT_THAT(xDer.XZ(), DoubleEq(0));
  ASSERT_THAT(xDer.YZ(), DoubleEq(0));
}

TEST_F(AVectorDerivatives3DTest, IsCorrectForNorm) {
  AutomaticDifferentiation::Second3D X(5 * x, 5, 0, 0);
  AutomaticDifferentiation::Second3D Y(2 * y, 0, 2, 0);
  AutomaticDifferentiation::Second3D Z(-z, 0, 0, -1);

  VectorDerivatives3D dv(X, Y, Z);

  auto norm = dv.norm();
  double root = std::sqrt(25 * x * x + 4 * y * y + z * z);

  ASSERT_THAT(norm.value(), DoubleEq(root));
  ASSERT_THAT(norm.dx(), DoubleEq(5 * (5 * x) / root));
  ASSERT_THAT(norm.dy(), DoubleEq(2 * (2 * y) / root));
  ASSERT_THAT(norm.dz(), DoubleEq((-1) * (-z) / root));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */