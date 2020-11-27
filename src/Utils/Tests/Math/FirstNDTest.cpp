/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Math/AutomaticDifferentiation/FirstND.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
using namespace AutomaticDifferentiation;
namespace Tests {

// Test class AFirstNDClass
// This class tests the FirstND class, which handles first derivatives in N dimensions.
class AFirstNDClass : public Test {
 public:
  // Declare a FirstND object for the use in the tests.
  FirstND d1;

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

// Tests that the FirstND object is of zero-th dimension by default.
TEST_F(AFirstNDClass, IsOfDegreeZeroByDefault) {
  ASSERT_THAT(d1.value(), DoubleEq(0.0));
  ASSERT_THAT(d1.dimensions(), Eq(0));
}

// Tests the constructor that takes the value as a double and the derivative as an Eigen::Vector3d.
TEST_F(AFirstNDClass, CanBeInitializedWithValueAndVector) {
  FirstND d2(arbitraryV, Eigen::Vector3d(arbitraryD1, arbitraryD2, arbitraryD3));

  ASSERT_THAT(d2.value(), DoubleEq(arbitraryV));
  ASSERT_THAT(d2.derivative(0), DoubleEq(arbitraryD1));
  ASSERT_THAT(d2.derivative(1), DoubleEq(arbitraryD2));
  ASSERT_THAT(d2.derivative(2), DoubleEq(arbitraryD3));
}

// Tests that a minus sign in front of a FirstND object inverts the sign of the value and all derivatives.
TEST_F(AFirstNDClass, UnaryMinusWorks) {
  FirstND d2(arbitraryV, Eigen::Vector3d(arbitraryD1, arbitraryD2, arbitraryD3));

  FirstND d4 = -d2;

  ASSERT_THAT(d4.value(), DoubleEq(-arbitraryV));
  ASSERT_THAT(d4.derivative(0), DoubleEq(-arbitraryD1));
  ASSERT_THAT(d4.derivative(1), DoubleEq(-arbitraryD2));
  ASSERT_THAT(d4.derivative(2), DoubleEq(-arbitraryD3));
}

// Tests that two FirstND objects can be added.
TEST_F(AFirstNDClass, CanBeAddedWithAnotherOne) {
  FirstND d2(arbitraryV, Eigen::Vector3d(arbitraryD1, arbitraryD2, arbitraryD3));
  FirstND d3(arbitraryVV, Eigen::Vector3d(arbitraryDD1, arbitraryDD2, arbitraryDD3));

  FirstND d4 = d2 + d3;
  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV + arbitraryVV));
  ASSERT_THAT(d4.derivative(0), DoubleEq(arbitraryD1 + arbitraryDD1));
  ASSERT_THAT(d4.derivative(1), DoubleEq(arbitraryD2 + arbitraryDD2));
  ASSERT_THAT(d4.derivative(2), DoubleEq(arbitraryD3 + arbitraryDD3));
}

// Tests, that a FirstND object can be multiplied by a scalar.
TEST_F(AFirstNDClass, CanBeMultipliedWithAScalar) {
  FirstND d2(arbitraryV, Eigen::Vector3d(arbitraryD1, arbitraryD2, arbitraryD3));
  double scalar = -2.04969;

  // Multiply from the left and from the right
  FirstND d4 = d2 * scalar;
  FirstND d5 = scalar * d2;

  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV * scalar));
  ASSERT_THAT(d4.derivative(0), DoubleEq(arbitraryD1 * scalar));
  ASSERT_THAT(d4.derivative(1), DoubleEq(arbitraryD2 * scalar));
  ASSERT_THAT(d4.derivative(2), DoubleEq(arbitraryD3 * scalar));
  ASSERT_THAT(d5.value(), DoubleEq(arbitraryV * scalar));
  ASSERT_THAT(d5.derivative(0), DoubleEq(arbitraryD1 * scalar));
  ASSERT_THAT(d5.derivative(1), DoubleEq(arbitraryD2 * scalar));
  ASSERT_THAT(d5.derivative(2), DoubleEq(arbitraryD3 * scalar));
}

// Tests that a FirstND object can be divided by a scalar.
TEST_F(AFirstNDClass, CanDivideByScalar) {
  FirstND d2(arbitraryV, Eigen::Vector3d(arbitraryD1, arbitraryD2, arbitraryD3));
  double scalar = -2.04969;

  FirstND d4 = d2 / scalar;

  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV / scalar));
  ASSERT_THAT(d4.derivative(0), DoubleEq(arbitraryD1 / scalar));
  ASSERT_THAT(d4.derivative(1), DoubleEq(arbitraryD2 / scalar));
  ASSERT_THAT(d4.derivative(2), DoubleEq(arbitraryD3 / scalar));
}

// Tests that two FirstND objects can be multiplied using the product rule.
TEST_F(AFirstNDClass, CanMultiplyWithAnotherDeriv) {
  FirstND d2(arbitraryV, Eigen::Vector3d(arbitraryD1, arbitraryD2, arbitraryD3));
  FirstND d3(arbitraryVV, Eigen::Vector3d(arbitraryDD1, arbitraryDD2, arbitraryDD3));

  FirstND d4 = d2 * d3;

  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV * arbitraryVV));
  ASSERT_THAT(d4.derivative(0), DoubleEq(arbitraryV * arbitraryDD1 + arbitraryD1 * arbitraryVV));
  ASSERT_THAT(d4.derivative(1), DoubleEq(arbitraryV * arbitraryDD2 + arbitraryD2 * arbitraryVV));
  ASSERT_THAT(d4.derivative(2), DoubleEq(arbitraryV * arbitraryDD3 + arbitraryD3 * arbitraryVV));
}

// Tests that two FirstND objects can be divided using the quotient rule.
TEST_F(AFirstNDClass, CanBeDividedByAnotherDeriv) {
  FirstND d2(arbitraryV, Eigen::Vector3d(arbitraryD1, arbitraryD2, arbitraryD3));
  FirstND d3(arbitraryVV, Eigen::Vector3d(arbitraryDD1, arbitraryDD2, arbitraryDD3));

  FirstND d4 = d2 / d3;

  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV / arbitraryVV));
  ASSERT_THAT(d4.derivative(0), DoubleEq(arbitraryD1 / arbitraryVV - arbitraryV / (arbitraryVV * arbitraryVV) * arbitraryDD1));
  ASSERT_THAT(d4.derivative(1), DoubleEq(arbitraryD2 / arbitraryVV - arbitraryV / (arbitraryVV * arbitraryVV) * arbitraryDD2));
  ASSERT_THAT(d4.derivative(2), DoubleEq(arbitraryD3 / arbitraryVV - arbitraryV / (arbitraryVV * arbitraryVV) * arbitraryDD3));
}

// Tests that a multiplication of two FirstND objects leads to the correct result.
TEST_F(AFirstNDClass, CorrectForProductXZ) {
  double x = 1.2568;
  double z = -4.2157;
  FirstND X(x, Eigen::Vector3d::UnitX());
  FirstND Z(z, Eigen::Vector3d::UnitZ());

  FirstND r = X * Z;

  ASSERT_THAT(r.value(), DoubleEq(x * z));
  ASSERT_THAT(r.derivative(0), DoubleEq(z));
  ASSERT_THAT(r.derivative(1), DoubleEq(0));
  ASSERT_THAT(r.derivative(2), DoubleEq(x));
}

// Tests that a multiplication of two First3D objects leads to the correct result.
TEST_F(AFirstNDClass, CorrectForProductX2_xyz) {
  double x = 1.2568;
  double y = 2.3131;
  double z = -4.2157;
  FirstND X2(x * x, Eigen::Vector3d(2 * x, 0, 0));
  FirstND XYZ(x * y * z, Eigen::Vector3d(y * z, x * z, x * y));

  FirstND r = X2 * XYZ;

  ASSERT_THAT(r.value(), DoubleEq(x * x * x * y * z));
  ASSERT_THAT(r.derivative(0), DoubleEq(3 * x * x * y * z));
  ASSERT_THAT(r.derivative(1), DoubleEq(x * x * x * z));
  ASSERT_THAT(r.derivative(2), DoubleEq(x * x * x * y));
}

// Tests that the correct result is given for a division of a FirstND object by itself.
TEST_F(AFirstNDClass, CorrectForDivision_X_X) {
  double x = 1.2568;
  FirstND X(x, Eigen::Vector3d::UnitX());

  FirstND r = X / X;

  ASSERT_THAT(r.value(), DoubleEq(1));
  ASSERT_THAT(r.derivative(0), DoubleEq(0));
  ASSERT_THAT(r.derivative(1), DoubleEq(0));
  ASSERT_THAT(r.derivative(2), DoubleEq(0));
}

// Tests, that the division of two FirstND objects leads to the correct result.
TEST_F(AFirstNDClass, CorrectForDivision_X_Z) {
  double x = 1.2568;
  double z = -4.2157;
  FirstND X(x, Eigen::Vector3d::UnitX());
  FirstND Z(z, Eigen::Vector3d::UnitZ());

  FirstND r = X / Z;

  ASSERT_THAT(r.value(), DoubleEq(x / z));
  ASSERT_THAT(r.derivative(0), DoubleEq(1.0 / z));
  ASSERT_THAT(r.derivative(1), DoubleEq(0));
  ASSERT_THAT(r.derivative(2), DoubleEq(-x / (z * z)));
}

// Tests, that the division of two FirstND objects leads to the correct result.
TEST_F(AFirstNDClass, CorrectForDivision_X2_xyz) {
  double x = 1.2568;
  double y = 2.3131;
  double z = -4.2157;
  FirstND X2(x * x, Eigen::Vector3d(2 * x, 0, 0));
  FirstND XYZ(x * y * z, Eigen::Vector3d(y * z, x * z, x * y));

  FirstND r = X2 / XYZ;

  ASSERT_THAT(r.value(), DoubleEq(x / (y * z)));
  ASSERT_THAT(r.derivative(0), DoubleEq(1 / (y * z)));
  ASSERT_THAT(r.derivative(1), DoubleEq(-x / (y * y * z)));
  ASSERT_THAT(r.derivative(2), DoubleEq(-x / (y * z * z)));
}

// Tests, that the square root function is correctly implemented for the FirstND object.
TEST_F(AFirstNDClass, CorrectForSqrt_x2y2z2) {
  double x = 1.2568;
  double y = 2.3131;
  double z = -4.2157;
  FirstND X(x, Eigen::Vector3d::UnitX());
  FirstND Y(y, Eigen::Vector3d::UnitY());
  FirstND Z(z, Eigen::Vector3d::UnitZ());

  FirstND R = sqrt(X * X + Y * Y + Z * Z);

  double r = std::sqrt(x * x + y * y + z * z);

  ASSERT_THAT(R.value(), DoubleEq(r));
  ASSERT_THAT(R.derivative(0), DoubleEq(x / r));
  ASSERT_THAT(R.derivative(1), DoubleEq(y / r));
  ASSERT_THAT(R.derivative(2), DoubleEq(z / r));
}

// Tests, that the exponential function is correctly implemented for the FirstND object.
TEST_F(AFirstNDClass, CorrectForExp_x2y) {
  double x = 1.2568;
  double y = 2.3131;
  FirstND X(x, Eigen::Vector3d::UnitX());
  FirstND Y(y, Eigen::Vector3d::UnitY());

  FirstND X2Y = X * X * Y;

  FirstND E = exp(X2Y);
  double expV = std::exp(x * x * y);

  ASSERT_THAT(E.value(), DoubleEq(expV));
  ASSERT_THAT(E.derivative(0), DoubleEq(expV * 2 * x * y));
  ASSERT_THAT(E.derivative(1), DoubleEq(expV * x * x));
  ASSERT_THAT(E.derivative(2), DoubleEq(0));
}

// Tests, that the square function and the exponential function are correctly implemented for the FirstND object.
TEST_F(AFirstNDClass, CorrectForSquare_Exp_x2y) {
  double x = 1.2568;
  double y = 2.3131;
  FirstND X(x, Eigen::Vector3d::UnitX());
  FirstND Y(y, Eigen::Vector3d::UnitY());

  FirstND X2Y = X * X * Y;

  // First use the exponential function.
  FirstND E = exp(X2Y);
  // Use the square function on the exponential function.
  FirstND Square = square(E);
  double expV = std::exp(x * x * y);
  double expV2 = expV * expV;

  ASSERT_THAT(Square.value(), DoubleEq(expV2));
  ASSERT_THAT(Square.derivative(0), DoubleEq(2 * expV2 * 2 * x * y));
  ASSERT_THAT(Square.derivative(1), DoubleEq(2 * expV2 * x * x));
  ASSERT_THAT(Square.derivative(2), DoubleEq(0));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
