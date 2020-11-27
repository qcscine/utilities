/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
using namespace AutomaticDifferentiation;
namespace Tests {

// Test class ADerivClass
// This class tests the First3D class, which handles first derivatives in 3 dimensions.
class ADerivClass : public Test {
 public:
  // Declare a First3D object for the use in the tests.
  First3D d1;

  // Set some arbitrary values for some function values and derivatives, which are used in the tests.
  double arbitraryV, arbitraryD1, arbitraryD2, arbitraryD3;
  double arbitraryVV, arbitraryDD1, arbitraryDD2, arbitraryDD3;

  void SetUp() override {
    arbitraryV = 0.23;
    arbitraryD1 = 1.111;
    arbitraryD2 = -9.23;
    arbitraryD3 = 29.3;
    arbitraryVV = 8.223;
    arbitraryDD1 = -2.111;
    arbitraryDD2 = 3.23;
    arbitraryDD3 = 9.0003;
  }
};

// Test that the First3D object d1 has its value and derivatives equal to zero by default.
TEST_F(ADerivClass, ContainsZeroElementsByDefault) {
  ASSERT_THAT(d1.value(), DoubleEq(0.0));
  ASSERT_THAT(d1.derivatives().norm(), DoubleEq(0.0));
}

// Tests the constructor of First3D that takes the value and the derivatives as doubles.
TEST_F(ADerivClass, CanBeInitializedWithValues) {
  First3D d2(arbitraryV, arbitraryD1, arbitraryD2, arbitraryD3);

  ASSERT_THAT(d2.value(), DoubleEq(arbitraryV));
  ASSERT_THAT(d2.derivatives()(0), DoubleEq(arbitraryD1));
  ASSERT_THAT(d2.derivatives()(1), DoubleEq(arbitraryD2));
  ASSERT_THAT(d2.derivatives()(2), DoubleEq(arbitraryD3));
}

// Tests the constructor of First3D that takes the value as a double and the derivatives as an Eigen::Vector3d.
TEST_F(ADerivClass, CanBeInitializedWithValueAndVector) {
  First3D d2(arbitraryV, Eigen::Vector3d(arbitraryD1, arbitraryD2, arbitraryD3));

  ASSERT_THAT(d2.value(), DoubleEq(arbitraryV));
  ASSERT_THAT(d2.derivatives()(0), DoubleEq(arbitraryD1));
  ASSERT_THAT(d2.derivatives()(1), DoubleEq(arbitraryD2));
  ASSERT_THAT(d2.derivatives()(2), DoubleEq(arbitraryD3));
}

// Tests that a minus sign in front of a First3D object inverts the sign of the value and all derivatives.
TEST_F(ADerivClass, UnaryMinusWorks) {
  First3D d2(arbitraryV, Eigen::Vector3d(arbitraryD1, arbitraryD2, arbitraryD3));

  First3D d4 = -d2;

  ASSERT_THAT(d4.value(), DoubleEq(-arbitraryV));
  ASSERT_THAT(d4.derivatives()(0), DoubleEq(-arbitraryD1));
  ASSERT_THAT(d4.derivatives()(1), DoubleEq(-arbitraryD2));
  ASSERT_THAT(d4.derivatives()(2), DoubleEq(-arbitraryD3));
}

// Tests that two First3D objects can be added.
TEST_F(ADerivClass, CanBeAddedWithAnotherOne) {
  First3D d2(arbitraryV, arbitraryD1, arbitraryD2, arbitraryD3);
  First3D d3(arbitraryVV, arbitraryDD1, arbitraryDD2, arbitraryDD3);

  First3D d4 = d2 + d3;
  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV + arbitraryVV));
  ASSERT_THAT(d4.derivatives()(0), DoubleEq(arbitraryD1 + arbitraryDD1));
  ASSERT_THAT(d4.derivatives()(1), DoubleEq(arbitraryD2 + arbitraryDD2));
  ASSERT_THAT(d4.derivatives()(2), DoubleEq(arbitraryD3 + arbitraryDD3));
}

// Tests, that a First3D object can be multiplied by a scalar.
TEST_F(ADerivClass, CanBeMultipliedWithAScalar) {
  First3D d2(arbitraryV, arbitraryD1, arbitraryD2, arbitraryD3);
  double scalar = -2.04969;

  // multiply from the left and the right
  First3D d4 = d2 * scalar;
  First3D d5 = scalar * d2;

  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV * scalar));
  ASSERT_THAT(d4.derivatives()(0), DoubleEq(arbitraryD1 * scalar));
  ASSERT_THAT(d4.derivatives()(1), DoubleEq(arbitraryD2 * scalar));
  ASSERT_THAT(d4.derivatives()(2), DoubleEq(arbitraryD3 * scalar));
  ASSERT_THAT(d5.value(), DoubleEq(arbitraryV * scalar));
  ASSERT_THAT(d5.derivatives()(0), DoubleEq(arbitraryD1 * scalar));
  ASSERT_THAT(d5.derivatives()(1), DoubleEq(arbitraryD2 * scalar));
  ASSERT_THAT(d5.derivatives()(2), DoubleEq(arbitraryD3 * scalar));
}

// Tests that a First3D object can be divided by a scalar.
TEST_F(ADerivClass, CanDivideByScalar) {
  First3D d2(arbitraryV, arbitraryD1, arbitraryD2, arbitraryD3);
  double scalar = -2.04969;

  First3D d4 = d2 / scalar;

  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV / scalar));
  ASSERT_THAT(d4.derivatives()(0), DoubleEq(arbitraryD1 / scalar));
  ASSERT_THAT(d4.derivatives()(1), DoubleEq(arbitraryD2 / scalar));
  ASSERT_THAT(d4.derivatives()(2), DoubleEq(arbitraryD3 / scalar));
}

// Tests that two First3D objects can be multiplied using the product rule.
TEST_F(ADerivClass, CanMultiplyWithAnotherDeriv) {
  First3D d2(arbitraryV, arbitraryD1, arbitraryD2, arbitraryD3);
  First3D d3(arbitraryVV, arbitraryDD1, arbitraryDD2, arbitraryDD3);

  First3D d4 = d2 * d3;

  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV * arbitraryVV));
  ASSERT_THAT(d4.derivatives()(0), DoubleEq(arbitraryV * arbitraryDD1 + arbitraryD1 * arbitraryVV));
  ASSERT_THAT(d4.derivatives()(1), DoubleEq(arbitraryV * arbitraryDD2 + arbitraryD2 * arbitraryVV));
  ASSERT_THAT(d4.derivatives()(2), DoubleEq(arbitraryV * arbitraryDD3 + arbitraryD3 * arbitraryVV));
}

// Tests that two First3D objects can be divided using the quotient rule.
TEST_F(ADerivClass, CanBeDividedByAnotherDeriv) {
  First3D d2(arbitraryV, arbitraryD1, arbitraryD2, arbitraryD3);
  First3D d3(arbitraryVV, arbitraryDD1, arbitraryDD2, arbitraryDD3);

  First3D d4 = d2 / d3;

  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV / arbitraryVV));
  ASSERT_THAT(d4.derivatives()(0),
              DoubleEq(arbitraryD1 / arbitraryVV - arbitraryV / (arbitraryVV * arbitraryVV) * arbitraryDD1));
  ASSERT_THAT(d4.derivatives()(1),
              DoubleEq(arbitraryD2 / arbitraryVV - arbitraryV / (arbitraryVV * arbitraryVV) * arbitraryDD2));
  ASSERT_THAT(d4.derivatives()(2),
              DoubleEq(arbitraryD3 / arbitraryVV - arbitraryV / (arbitraryVV * arbitraryVV) * arbitraryDD3));
}

// Tests that a multiplication of two First3D objects leads to the correct result.
TEST_F(ADerivClass, CorrectForProductXZ) {
  double x = 1.2568;
  double z = -4.2157;
  First3D X(x, 1, 0, 0);
  First3D Z(z, 0, 0, 1);

  First3D r = X * Z;

  // The value should be just the product of the doubles x and z.
  ASSERT_THAT(r.value(), DoubleEq(x * z));
  // The resulting derivative should be (z, 0, x).
  ASSERT_THAT(r.derivatives()(0), DoubleEq(z));
  ASSERT_THAT(r.derivatives()(1), DoubleEq(0));
  ASSERT_THAT(r.derivatives()(2), DoubleEq(x));
}

// Tests that a multiplication of two First3D objects leads to the correct result.
TEST_F(ADerivClass, CorrectForProductX2_xyz) {
  double x = 1.2568;
  double y = 2.3131;
  double z = -4.2157;
  First3D X2(x * x, 2 * x, 0, 0);
  First3D XYZ(x * y * z, y * z, x * z, x * y);

  First3D r = X2 * XYZ;

  ASSERT_THAT(r.value(), DoubleEq(x * x * x * y * z));
  ASSERT_THAT(r.derivatives()(0), DoubleEq(3 * x * x * y * z));
  ASSERT_THAT(r.derivatives()(1), DoubleEq(x * x * x * z));
  ASSERT_THAT(r.derivatives()(2), DoubleEq(x * x * x * y));
}

// Tests that the correct result is given for a division of a First3D object by itself.
TEST_F(ADerivClass, CorrectForDivision_X_X) {
  double x = 1.2568;
  First3D X(x, 1, 0, 0);

  First3D r = X / X;

  ASSERT_THAT(r.value(), DoubleEq(1));
  ASSERT_THAT(r.derivatives()(0), DoubleEq(0));
  ASSERT_THAT(r.derivatives()(1), DoubleEq(0));
  ASSERT_THAT(r.derivatives()(2), DoubleEq(0));
}

// Tests, that the division of two First3D objects leads to the correct result.
TEST_F(ADerivClass, CorrectForDivision_X_Z) {
  double x = 1.2568;
  double z = -4.2157;
  First3D X(x, 1, 0, 0);
  First3D Z(z, 0, 0, 1);

  First3D r = X / Z;

  ASSERT_THAT(r.value(), DoubleEq(x / z));
  ASSERT_THAT(r.derivatives()(0), DoubleEq(1.0 / z));
  ASSERT_THAT(r.derivatives()(1), DoubleEq(0));
  ASSERT_THAT(r.derivatives()(2), DoubleEq(-x / (z * z)));
}

// Tests, that the division of two First3D objects leads to the correct result.
TEST_F(ADerivClass, CorrectForDivision_X2_xyz) {
  double x = 1.2568;
  double y = 2.3131;
  double z = -4.2157;
  First3D X2(x * x, 2 * x, 0, 0);
  First3D XYZ(x * y * z, y * z, x * z, x * y);

  First3D r = X2 / XYZ;

  ASSERT_THAT(r.value(), DoubleEq(x / (y * z)));
  ASSERT_THAT(r.derivatives()(0), DoubleEq(1 / (y * z)));
  ASSERT_THAT(r.derivatives()(1), DoubleEq(-x / (y * y * z)));
  ASSERT_THAT(r.derivatives()(2), DoubleEq(-x / (y * z * z)));
}

// Tests, that the square root function is correctly implemented for the First3D object.
TEST_F(ADerivClass, CorrectForSqrt_x2y2z2) {
  double x = 1.2568;
  double y = 2.3131;
  double z = -4.2157;
  First3D X(x, 1, 0, 0);
  First3D Y(y, 0, 1, 0);
  First3D Z(z, 0, 0, 1);

  First3D R = sqrt(X * X + Y * Y + Z * Z);

  double r = std::sqrt(x * x + y * y + z * z);

  ASSERT_THAT(R.value(), DoubleEq(r));
  ASSERT_THAT(R.derivatives()(0), DoubleEq(x / r));
  ASSERT_THAT(R.derivatives()(1), DoubleEq(y / r));
  ASSERT_THAT(R.derivatives()(2), DoubleEq(z / r));
}

// Tests, that the exponential function is correctly implemented for the First3D object.
TEST_F(ADerivClass, CorrectForExp_x2y) {
  double x = 1.2568;
  double y = 2.3131;
  First3D X(x, 1, 0, 0);
  First3D Y(y, 0, 1, 0);

  First3D X2Y = X * X * Y;

  First3D E = exp(X2Y);
  double expV = std::exp(x * x * y);

  ASSERT_THAT(E.value(), DoubleEq(expV));
  ASSERT_THAT(E.derivatives()(0), DoubleEq(expV * 2 * x * y));
  ASSERT_THAT(E.derivatives()(1), DoubleEq(expV * x * x));
  ASSERT_THAT(E.derivatives()(2), DoubleEq(0));
}

// Tests that the function, that creates a First3D object from the one-dimensional First1D object,
// gives the correct result.
TEST_F(ADerivClass, CanCreate3DFrom1DInXDirection) {
  double v = 0.4432, d = -2.213;
  double x = 9.232;
  // The direction of the one-dimensional object is along the x-axis.
  Eigen::Vector3d vec(x, 0, 0);

  First1D d1D(v, d);
  auto d3D = get3Dfrom1D<DerivativeOrder::One>(d1D, vec);

  const auto& deriv = d3D.derivatives();
  ASSERT_THAT(d3D.value(), DoubleEq(d1D.value()));
  ASSERT_THAT(deriv.x(), DoubleEq(d1D.derivatives()));
  ASSERT_THAT(deriv.y(), DoubleEq(0));
  ASSERT_THAT(deriv.z(), DoubleEq(0));
}

// Tests that the function, that creates a First3D object from the one-dimensional First1D object,
// gives the correct result.
TEST_F(ADerivClass, CanCreate3DFrom1DInRandomDirection) {
  double v = 0.4432, d = -2.213;
  double x = 9.232, y = -0.23, z = 1.0;
  // The direction of the one-dimensional object is randomly chosen (not along the x, y or z-axis).
  Eigen::Vector3d vec(x, y, z);
  double R = vec.norm();

  First1D d1D(v, d);
  auto d3D = get3Dfrom1D<DerivativeOrder::One>(d1D, vec);

  const auto& deriv = d3D.derivatives();
  ASSERT_THAT(d3D.value(), DoubleEq(d1D.value()));
  ASSERT_THAT(deriv.x(), DoubleEq(d1D.derivatives() * x / R));
  ASSERT_THAT(deriv.y(), DoubleEq(d1D.derivatives() * y / R));
  ASSERT_THAT(deriv.z(), DoubleEq(d1D.derivatives() * z / R));
}

// Tests that a double can be divided by a First3D object correctly.
TEST_F(ADerivClass, DoubleCanBeDividedByValueWithDerivatives) {
  double v = 0.4432, d = -2.213;
  double x = -9.3343;

  First1D d1D(v, d);
  auto xWithDerivatives = constant1D<DerivativeOrder::One>(x);

  auto expected = xWithDerivatives / d1D;
  auto obtained = x / d1D;

  ASSERT_THAT(obtained.value(), DoubleEq(expected.value()));
  ASSERT_THAT(obtained.derivatives(), DoubleNear(expected.derivatives(), 1e-10));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
