/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Utils/Math/AutomaticDifferentiation/Second3D.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
using namespace AutomaticDifferentiation;
namespace Tests {

// Test class ASecond3DClass
// This class tests the Second3D class, which handles first and second derivatives in 3 dimensions.
class ASecond3DClass : public Test {
 public:
  // Declare a Second3D object for the use in the tests.
  Second3D d1;

  // Set some arbitrary values for some function values and derivatives, which are used in the tests.
  double arbitraryV, arbitraryD1, arbitraryD2, arbitraryD3, arbitraryXX, arbitraryYY, arbitraryZZ, arbitraryXY,
      arbitraryXZ, arbitraryYZ;
  double arbitraryVV, arbitraryDD1, arbitraryDD2, arbitraryDD3, arbitraryXX2, arbitraryYY2, arbitraryZZ2, arbitraryXY2,
      arbitraryXZ2, arbitraryYZ2;

  void SetUp() override {
    arbitraryV = 0.23;
    arbitraryD1 = 1.111;
    arbitraryD2 = -9.23;
    arbitraryD3 = 29.3;
    arbitraryVV = 8.223;
    arbitraryDD1 = -2.111;
    arbitraryDD2 = 3.23;
    arbitraryDD3 = 9.0003;
    arbitraryXX = 2.231;
    arbitraryYY = 0.22;
    arbitraryZZ = -5.02;
    arbitraryXY = -0.1234;
    arbitraryXZ = 1.2222;
    arbitraryYZ = 3.025;
    arbitraryXX2 = -0.0125;
    arbitraryYY2 = -5.45;
    arbitraryZZ2 = 2.222;
    arbitraryXY2 = 1.0001;
    arbitraryXZ2 = -0.9857;
    arbitraryYZ2 = 0.665;
  }
};

// Tests that the value, the first derivatives and second derivatives are all zero by default.
TEST_F(ASecond3DClass, ContainsZeroElementsByDefault) {
  ASSERT_THAT(d1.value(), DoubleEq(0.0));
  ASSERT_THAT(d1.dx(), DoubleEq(0.0));
  ASSERT_THAT(d1.dy(), DoubleEq(0.0));
  ASSERT_THAT(d1.dz(), DoubleEq(0.0));
  ASSERT_THAT(d1.XX(), DoubleEq(0.0));
  ASSERT_THAT(d1.YY(), DoubleEq(0.0));
  ASSERT_THAT(d1.ZZ(), DoubleEq(0.0));
  ASSERT_THAT(d1.XY(), DoubleEq(0.0));
  ASSERT_THAT(d1.XZ(), DoubleEq(0.0));
  ASSERT_THAT(d1.YZ(), DoubleEq(0.0));
}

// Tests the constructor that takes the function value, the first derivatives and the second derivatives as doubles.
TEST_F(ASecond3DClass, CanBeInitializedWithValues) {
  Second3D d2(arbitraryV, arbitraryD1, arbitraryD2, arbitraryD3, arbitraryXX, arbitraryYY, arbitraryZZ, arbitraryXY,
              arbitraryXZ, arbitraryYZ);

  ASSERT_THAT(d2.value(), DoubleEq(arbitraryV));
  ASSERT_THAT(d2.dx(), DoubleEq(arbitraryD1));
  ASSERT_THAT(d2.dy(), DoubleEq(arbitraryD2));
  ASSERT_THAT(d2.dz(), DoubleEq(arbitraryD3));
  ASSERT_THAT(d2.XX(), DoubleEq(arbitraryXX));
  ASSERT_THAT(d2.YY(), DoubleEq(arbitraryYY));
  ASSERT_THAT(d2.ZZ(), DoubleEq(arbitraryZZ));
  ASSERT_THAT(d2.XY(), DoubleEq(arbitraryXY));
  ASSERT_THAT(d2.XZ(), DoubleEq(arbitraryXZ));
  ASSERT_THAT(d2.YZ(), DoubleEq(arbitraryYZ));
}

// Tests that a minus sign in front of a Second3D object changes the sign of all its components.
TEST_F(ASecond3DClass, UnaryMinusWorks) {
  Second3D d2(arbitraryV, arbitraryD1, arbitraryD2, arbitraryD3, arbitraryXX, arbitraryYY, arbitraryZZ, arbitraryXY,
              arbitraryXZ, arbitraryYZ);

  Second3D d4 = -d2;

  ASSERT_THAT(d4.value(), DoubleEq(-arbitraryV));
  ASSERT_THAT(d4.dx(), DoubleEq(-arbitraryD1));
  ASSERT_THAT(d4.dy(), DoubleEq(-arbitraryD2));
  ASSERT_THAT(d4.dz(), DoubleEq(-arbitraryD3));
  ASSERT_THAT(d4.XX(), DoubleEq(-arbitraryXX));
  ASSERT_THAT(d4.YY(), DoubleEq(-arbitraryYY));
  ASSERT_THAT(d4.ZZ(), DoubleEq(-arbitraryZZ));
  ASSERT_THAT(d4.XY(), DoubleEq(-arbitraryXY));
  ASSERT_THAT(d4.XZ(), DoubleEq(-arbitraryXZ));
  ASSERT_THAT(d4.YZ(), DoubleEq(-arbitraryYZ));
}

// Tests that two Second3D objects can be added.
TEST_F(ASecond3DClass, CanBeAddedWithAnotherOne) {
  Second3D d2(arbitraryV, arbitraryD1, arbitraryD2, arbitraryD3, arbitraryXX, arbitraryYY, arbitraryZZ, arbitraryXY,
              arbitraryXZ, arbitraryYZ);
  Second3D d3(arbitraryVV, arbitraryDD1, arbitraryDD2, arbitraryDD3, arbitraryXX2, arbitraryYY2, arbitraryZZ2,
              arbitraryXY2, arbitraryXZ2, arbitraryYZ2);

  Second3D d4 = d2 + d3;

  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV + arbitraryVV));
  ASSERT_THAT(d4.dx(), DoubleEq(arbitraryD1 + arbitraryDD1));
  ASSERT_THAT(d4.dy(), DoubleEq(arbitraryD2 + arbitraryDD2));
  ASSERT_THAT(d4.dz(), DoubleEq(arbitraryD3 + arbitraryDD3));
  ASSERT_THAT(d4.XX(), DoubleEq(arbitraryXX + arbitraryXX2));
  ASSERT_THAT(d4.YY(), DoubleEq(arbitraryYY + arbitraryYY2));
  ASSERT_THAT(d4.ZZ(), DoubleEq(arbitraryZZ + arbitraryZZ2));
  ASSERT_THAT(d4.XY(), DoubleEq(arbitraryXY + arbitraryXY2));
  ASSERT_THAT(d4.XZ(), DoubleEq(arbitraryXZ + arbitraryXZ2));
  ASSERT_THAT(d4.YZ(), DoubleEq(arbitraryYZ + arbitraryYZ2));
}

// Tests that a Second3D object can be multiplied by a scalar.
TEST_F(ASecond3DClass, CanBeMultipliedWithAScalar) {
  Second3D d2(arbitraryV, arbitraryD1, arbitraryD2, arbitraryD3, arbitraryXX, arbitraryYY, arbitraryZZ, arbitraryXY,
              arbitraryXZ, arbitraryYZ);
  double scalar = -2.04969;

  // Multiply from the left and from the right
  Second3D d4 = d2 * scalar;
  Second3D d5 = scalar * d2;

  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV * scalar));
  ASSERT_THAT(d4.dx(), DoubleEq(arbitraryD1 * scalar));
  ASSERT_THAT(d4.dy(), DoubleEq(arbitraryD2 * scalar));
  ASSERT_THAT(d4.dz(), DoubleEq(arbitraryD3 * scalar));
  ASSERT_THAT(d4.XX(), DoubleEq(arbitraryXX * scalar));
  ASSERT_THAT(d4.YY(), DoubleEq(arbitraryYY * scalar));
  ASSERT_THAT(d4.ZZ(), DoubleEq(arbitraryZZ * scalar));
  ASSERT_THAT(d4.XY(), DoubleEq(arbitraryXY * scalar));
  ASSERT_THAT(d4.XZ(), DoubleEq(arbitraryXZ * scalar));
  ASSERT_THAT(d4.YZ(), DoubleEq(arbitraryYZ * scalar));

  ASSERT_THAT(d5.value(), DoubleEq(arbitraryV * scalar));
  ASSERT_THAT(d5.dx(), DoubleEq(arbitraryD1 * scalar));
  ASSERT_THAT(d5.dy(), DoubleEq(arbitraryD2 * scalar));
  ASSERT_THAT(d5.dz(), DoubleEq(arbitraryD3 * scalar));
  ASSERT_THAT(d5.XX(), DoubleEq(arbitraryXX * scalar));
  ASSERT_THAT(d5.YY(), DoubleEq(arbitraryYY * scalar));
  ASSERT_THAT(d5.ZZ(), DoubleEq(arbitraryZZ * scalar));
  ASSERT_THAT(d5.XY(), DoubleEq(arbitraryXY * scalar));
  ASSERT_THAT(d5.XZ(), DoubleEq(arbitraryXZ * scalar));
  ASSERT_THAT(d5.YZ(), DoubleEq(arbitraryYZ * scalar));
}

// Tests that a Second3D object can be divided by a scalar.
TEST_F(ASecond3DClass, CanDivideByScalar) {
  Second3D d2(arbitraryV, arbitraryD1, arbitraryD2, arbitraryD3, arbitraryXX, arbitraryYY, arbitraryZZ, arbitraryXY,
              arbitraryXZ, arbitraryYZ);
  double scalar = -2.04969;

  Second3D d4 = d2 / scalar;

  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV / scalar));
  ASSERT_THAT(d4.dx(), DoubleEq(arbitraryD1 / scalar));
  ASSERT_THAT(d4.dy(), DoubleEq(arbitraryD2 / scalar));
  ASSERT_THAT(d4.dz(), DoubleEq(arbitraryD3 / scalar));
  ASSERT_THAT(d4.XX(), DoubleEq(arbitraryXX / scalar));
  ASSERT_THAT(d4.YY(), DoubleEq(arbitraryYY / scalar));
  ASSERT_THAT(d4.ZZ(), DoubleEq(arbitraryZZ / scalar));
  ASSERT_THAT(d4.XY(), DoubleEq(arbitraryXY / scalar));
  ASSERT_THAT(d4.XZ(), DoubleEq(arbitraryXZ / scalar));
  ASSERT_THAT(d4.YZ(), DoubleEq(arbitraryYZ / scalar));
}

// Tests that a Second3D object is multiplied correctly by itself.
TEST_F(ASecond3DClass, CorrectForProductXX) {
  double x = 1.2568;
  Second3D X(x, 1, 0, 0);

  Second3D r = X * X;

  ASSERT_THAT(r.value(), DoubleEq(x * x));
  ASSERT_THAT(r.dx(), DoubleEq(2 * x));
  ASSERT_THAT(r.dy(), DoubleEq(0));
  ASSERT_THAT(r.dz(), DoubleEq(0));
  ASSERT_THAT(r.XX(), DoubleEq(2));
  ASSERT_THAT(r.YY(), DoubleEq(0));
  ASSERT_THAT(r.ZZ(), DoubleEq(0));
  ASSERT_THAT(r.XY(), DoubleEq(0));
  ASSERT_THAT(r.XZ(), DoubleEq(0));
  ASSERT_THAT(r.YZ(), DoubleEq(0));
}

// Tests that the functions x*y*z and x^2 are multiplied correctly.
TEST_F(ASecond3DClass, CorrectForProductXYZ_2X2) {
  double x = 1.2568;
  double y = -0.45;
  double z = 2.3456;
  Second3D XYZ(x * y * z, y * z, x * z, x * y, 0, 0, 0, z, y, x);
  Second3D X2(x * x, 2 * x, 0, 0, 2, 0, 0, 0, 0, 0);

  Second3D r = XYZ * X2;

  ASSERT_THAT(r.value(), DoubleEq(x * x * x * y * z));
  ASSERT_THAT(r.dx(), DoubleEq(3 * x * x * y * z));
  ASSERT_THAT(r.dy(), DoubleEq(x * x * x * z));
  ASSERT_THAT(r.dz(), DoubleEq(x * x * x * y));
  ASSERT_THAT(r.XX(), DoubleEq(6 * x * y * z));
  ASSERT_THAT(r.YY(), DoubleEq(0));
  ASSERT_THAT(r.ZZ(), DoubleEq(0));
  ASSERT_THAT(r.XY(), DoubleEq(3 * x * x * z));
  ASSERT_THAT(r.XZ(), DoubleEq(3 * x * x * y));
  ASSERT_THAT(r.YZ(), DoubleEq(x * x * x));
}

// Tests that the division of the function x*y*z by the function x^2 is performed correctly.
TEST_F(ASecond3DClass, CorrectForDivision_XYZ_X2) {
  double x = 1.2568;
  double y = -0.45;
  double z = 2.3456;
  Second3D XYZ(x * y * z, y * z, x * z, x * y, 0, 0, 0, z, y, x);
  Second3D X2(x * x, 2 * x, 0, 0, 2, 0, 0, 0, 0, 0);

  Second3D r = XYZ / X2;

  ASSERT_THAT(r.value(), DoubleEq(y * z / x));
  ASSERT_THAT(r.dx(), DoubleEq(-y * z / (x * x)));
  ASSERT_THAT(r.dy(), DoubleEq(z / x));
  ASSERT_THAT(r.dz(), DoubleEq(y / x));
  ASSERT_THAT(r.XX(), DoubleEq(2 * y * z / (x * x * x)));
  ASSERT_THAT(r.YY(), DoubleEq(0));
  ASSERT_THAT(r.ZZ(), DoubleEq(0));
  ASSERT_THAT(r.XY(), DoubleEq(-z / (x * x)));
  ASSERT_THAT(r.XZ(), DoubleEq(-y / (x * x)));
  ASSERT_THAT(r.YZ(), DoubleEq(1.0 / x));
}

// Tests that the square root function is correctly implemented.
TEST_F(ASecond3DClass, CorrectForSqrt_x3y_divByZ) {
  double x = 1.2568;
  double y = 0.45;
  double z = 2.3456;
  Second3D X(x, 1, 0, 0, 0, 0, 0, 0, 0, 0);
  Second3D Y(y, 0, 1, 0, 0, 0, 0, 0, 0, 0);
  Second3D Z(z, 0, 0, 1, 0, 0, 0, 0, 0, 0);

  Second3D r = sqrt(X * X * X * Y / Z);

  double sqx = std::sqrt(x);
  double sqy = std::sqrt(y);
  double sqz = std::sqrt(z);

  ASSERT_THAT(r.value(), DoubleEq(sqx * sqx * sqx * sqy / sqz));
  ASSERT_THAT(r.dx(), DoubleEq(1.5 * sqx * sqy / sqz));
  ASSERT_THAT(r.dy(), DoubleEq(0.5 * sqx * sqx * sqx / (sqz * sqy)));
  ASSERT_THAT(r.dz(), DoubleEq(-0.5 * sqx * x * sqy / (z * sqz)));
  ASSERT_THAT(r.XX(), DoubleEq(1.5 * 0.5 * sqy / (sqz * sqx)));
  ASSERT_THAT(r.YY(), DoubleEq(-0.5 * 0.5 * x * sqx / (y * sqy * sqz)));
  ASSERT_THAT(r.ZZ(), DoubleEq(0.5 * 1.5 * x * sqx * sqy / (z * z * sqz)));
  ASSERT_THAT(r.XY(), DoubleNear(1.5 * 0.5 * sqx / (sqz * sqy), 1e-10));
  ASSERT_THAT(r.XZ(), DoubleNear(-1.5 * 0.5 * sqx * sqy / (z * sqz), 1e-10));
  ASSERT_THAT(r.YZ(), DoubleNear(-0.5 * 0.5 * x * sqx / (sqy * z * sqz), 1e-10));
}

// Tests that the exponential function is correctly implemented.
TEST_F(ASecond3DClass, CorrectForExp_xyz2ByX) {
  double x = 1.2568;
  double y = 0.45;
  double z = 2.3456;
  Second3D X(x, 1, 0, 0, 0, 0, 0, 0, 0, 0);
  Second3D Y(y, 0, 1, 0, 0, 0, 0, 0, 0, 0);
  Second3D Z(z, 0, 0, 1, 0, 0, 0, 0, 0, 0);

  Second3D r = exp(X * Y * Z * Z / X);

  double ex = std::exp(y * z * z);

  ASSERT_THAT(r.value(), DoubleEq(ex));
  ASSERT_THAT(r.dx(), DoubleNear(0, 1e-10));
  ASSERT_THAT(r.dy(), DoubleEq(ex * z * z));
  ASSERT_THAT(r.dz(), DoubleEq(ex * 2 * y * z));
  ASSERT_THAT(r.XX(), DoubleNear(0, 1e-10));
  ASSERT_THAT(r.YY(), DoubleEq(z * z * z * z * ex));
  ASSERT_THAT(r.ZZ(), DoubleEq(2 * y * (ex + z * ex * 2 * y * z)));
  ASSERT_THAT(r.XY(), DoubleNear(0, 1e-10));
  ASSERT_THAT(r.XZ(), DoubleNear(0, 1e-10));
  ASSERT_THAT(r.YZ(), DoubleEq(2 * z * ex + z * z * ex * y * 2 * z));
}

// Tests that the arccos function is correctly implemented.
TEST_F(ASecond3DClass, CorrectForArccos_xyz2ByX) {
  double x = -1.2568;
  double y = 0.45;
  double z = 1.3456;
  Second3D X(x, 1, 0, 0, 0, 0, 0, 0, 0, 0);
  Second3D Y(y, 0, 1, 0, 0, 0, 0, 0, 0, 0);
  Second3D Z(z, 0, 0, 1, 0, 0, 0, 0, 0, 0);

  Second3D r = arccos(X * Y * Z * Z / X);

  double ex = std::acos(y * z * z);
  double root = std::sqrt(1 - y * y * z * z * z * z);

  ASSERT_THAT(r.value(), DoubleEq(ex));
  ASSERT_THAT(r.dx(), DoubleNear(0, 1e-10));
  ASSERT_THAT(r.dy(), DoubleNear(-z * z / root, 1e-10));
  ASSERT_THAT(r.dz(), DoubleNear(-2 * y * z / root, 1e-10));
  ASSERT_THAT(r.XX(), DoubleNear(0, 1e-10));
  ASSERT_THAT(r.YY(), DoubleNear(-y * z * z * z * z * z * z / (root * root * root), 1e-10));
  ASSERT_THAT(r.ZZ(), DoubleNear(-4 * y * y * y * z * z * z * z / (root * root * root) - 2 * y / root, 1e-10));
  ASSERT_THAT(r.XY(), DoubleNear(0, 1e-10));
  ASSERT_THAT(r.XZ(), DoubleNear(0, 1e-10));
  ASSERT_THAT(r.YZ(), DoubleNear(-2 * y * y * z * z * z * z * z / (root * root * root) - 2 * z / root, 1e-10));
}

// Tests that a Second3D object can be obtained from a Second1D object correctly along the vector of the x-axis.
TEST_F(ASecond3DClass, CanCreate3DFrom1DInXDirection) {
  double v = 0.4432, d = -2.213, h = 1.34;
  double x = 9.232;
  Eigen::Vector3d vec(x, 0, 0);

  Second1D h1D(v, d, h);
  auto h3D = get3Dfrom1D<DerivativeOrder::Two>(h1D, vec);

  ASSERT_THAT(h3D.value(), DoubleEq(h1D.value()));
  ASSERT_THAT(h3D.dx(), DoubleEq(h1D.first()));
  ASSERT_THAT(h3D.dy(), DoubleEq(0));
  ASSERT_THAT(h3D.dz(), DoubleEq(0));
  ASSERT_THAT(h3D.XX(), DoubleEq(h1D.second()));
  ASSERT_THAT(h3D.YY(), DoubleEq(h1D.first() / x));
  ASSERT_THAT(h3D.ZZ(), DoubleEq(h1D.first() / x));
  ASSERT_THAT(h3D.XY(), DoubleEq(0));
  ASSERT_THAT(h3D.XZ(), DoubleEq(0));
  ASSERT_THAT(h3D.YZ(), DoubleEq(0));
}

// Tests that a Second3D object can be obtained from a Second1D object correctly along a random direction.
TEST_F(ASecond3DClass, CanCreate3DFrom1DInRandomDirection) {
  double v = 0.4432, d = -2.213, h = 1.34;
  double x = 9.232, y = -0.23, z = 1.0;
  double xx = x * x, yy = y * y, zz = z * z, xy = x * y, xz = x * z, yz = y * z;
  Eigen::Vector3d vec(x, y, z);
  double R = vec.norm();
  double RR = R * R, RRR = R * R * R;

  Second1D h1D(v, d, h);
  auto h3D = get3Dfrom1D<DerivativeOrder::Two>(h1D, vec);

  ASSERT_THAT(h3D.value(), DoubleEq(h1D.value()));
  ASSERT_THAT(h3D.dx(), DoubleEq(h1D.first() * x / R));
  ASSERT_THAT(h3D.dy(), DoubleEq(h1D.first() * y / R));
  ASSERT_THAT(h3D.dz(), DoubleEq(h1D.first() * z / R));
  ASSERT_THAT(h3D.XX(), DoubleEq(h1D.second() * xx / RR + h1D.first() * (1 / R - xx / RRR)));
  ASSERT_THAT(h3D.YY(), DoubleEq(h1D.second() * yy / RR + h1D.first() * (1 / R - yy / RRR)));
  ASSERT_THAT(h3D.ZZ(), DoubleEq(h1D.second() * zz / RR + h1D.first() * (1 / R - zz / RRR)));
  ASSERT_THAT(h3D.XY(), DoubleEq(h1D.second() * xy / RR + h1D.first() * (-xy / RRR)));
  ASSERT_THAT(h3D.XZ(), DoubleEq(h1D.second() * xz / RR + h1D.first() * (-xz / RRR)));
  ASSERT_THAT(h3D.YZ(), DoubleEq(h1D.second() * yz / RR + h1D.first() * (-yz / RRR)));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
