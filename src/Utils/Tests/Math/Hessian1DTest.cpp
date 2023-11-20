/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Math/AutomaticDifferentiation/Second1D.h>
#include <gmock/gmock.h>
#include <cmath>

using namespace testing;
namespace Scine {
namespace Utils {
using namespace AutomaticDifferentiation;
namespace Tests {

// Test class ASecond1DClass
// This class tests the Second1D class, which handles first and second derivatives in one dimension.
class ASecond1DClass : public Test {
 public:
  // Declare a Second1D object for the use in the tests.
  Second1D d1;

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

// Tests that the value, first derivative and second derivative is zero by default.
TEST_F(ASecond1DClass, ContainsZeroElementsByDefault) {
  ASSERT_THAT(d1.value(), DoubleEq(0.0));
  ASSERT_THAT(d1.first(), DoubleEq(0.0));
  ASSERT_THAT(d1.second(), DoubleEq(0.0));
}

// Tests the constructor that takes the function value, first derivative and second derivative arguments as doubles.
TEST_F(ASecond1DClass, CanBeInitializedWithValues) {
  Second1D d2(arbitraryV, arbitraryD1, arbitraryD2);

  ASSERT_THAT(d2.value(), DoubleEq(arbitraryV));
  ASSERT_THAT(d2.first(), DoubleEq(arbitraryD1));
  ASSERT_THAT(d2.second(), DoubleEq(arbitraryD2));
}

// Tests that a minus sign in front of a Second1D object changes the sign of all its components.
TEST_F(ASecond1DClass, UnaryMinusWorks) {
  Second1D d2(arbitraryV, arbitraryD1, arbitraryD2);

  Second1D d4 = -d2;

  ASSERT_THAT(d4.value(), DoubleEq(-arbitraryV));
  ASSERT_THAT(d4.first(), DoubleEq(-arbitraryD1));
  ASSERT_THAT(d4.second(), DoubleEq(-arbitraryD2));
}

// Tests that two Second1D objects can be added.
TEST_F(ASecond1DClass, CanBeAddedWithAnotherOne) {
  Second1D d2(arbitraryV, arbitraryD1, arbitraryD2);
  Second1D d3(arbitraryVV, arbitraryDD1, arbitraryDD2);

  Second1D d4 = d2 + d3;
  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV + arbitraryVV));
  ASSERT_THAT(d4.first(), DoubleEq(arbitraryD1 + arbitraryDD1));
  ASSERT_THAT(d4.second(), DoubleEq(arbitraryD2 + arbitraryDD2));
}

// Tests that a Second1D object can be multiplied by a scalar.
TEST_F(ASecond1DClass, CanBeMultipliedWithAScalar) {
  Second1D d2(arbitraryV, arbitraryD1, arbitraryD2);
  double scalar = -2.04969;

  // Multiply from the left and from the right
  Second1D d4 = d2 * scalar;
  Second1D d5 = scalar * d2;

  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV * scalar));
  ASSERT_THAT(d4.first(), DoubleEq(arbitraryD1 * scalar));
  ASSERT_THAT(d4.second(), DoubleEq(arbitraryD2 * scalar));
  ASSERT_THAT(d5.value(), DoubleEq(arbitraryV * scalar));
  ASSERT_THAT(d5.first(), DoubleEq(arbitraryD1 * scalar));
  ASSERT_THAT(d5.second(), DoubleEq(arbitraryD2 * scalar));
}

// Tests that a Second1D object can be divided by a scalar.
TEST_F(ASecond1DClass, CanDivideByScalar) {
  Second1D d2(arbitraryV, arbitraryD1, arbitraryD2);
  double scalar = -2.04969;

  Second1D d4 = d2 / scalar;

  ASSERT_THAT(d4.value(), DoubleEq(arbitraryV / scalar));
  ASSERT_THAT(d4.first(), DoubleEq(arbitraryD1 / scalar));
  ASSERT_THAT(d4.second(), DoubleEq(arbitraryD2 / scalar));
}

// Tests that a Second1D object is multiplied correctly by itself.
TEST_F(ASecond1DClass, CorrectForProductXX) {
  double x = 1.2568;
  Second1D X(x, 1, 0);

  Second1D r = X * X;

  ASSERT_THAT(r.value(), DoubleEq(x * x));
  ASSERT_THAT(r.first(), DoubleEq(2 * x));
  ASSERT_THAT(r.second(), DoubleEq(2));
}

// Tests that the function x^2 is correctly multiplied by the function 2*x.
TEST_F(ASecond1DClass, CorrectForProductX2_2X) {
  double x = 1.2568;
  Second1D X2(x * x, 2 * x, 2);
  Second1D X(2 * x, 2, 0);

  Second1D r = X2 * X;

  ASSERT_THAT(r.value(), DoubleEq(2 * x * x * x));
  ASSERT_THAT(r.first(), DoubleEq(6 * x * x));
  ASSERT_THAT(r.second(), DoubleEq(12 * x));
}

// Tests that a Second1D object is divided by itself correctly.
TEST_F(ASecond1DClass, CorrectForDivision_X_X) {
  double x = 1.2568;
  Second1D X(x, 1, 0);

  Second1D r = X / X;

  ASSERT_THAT(r.value(), DoubleEq(1));
  ASSERT_THAT(r.first(), DoubleEq(0));
  ASSERT_THAT(r.second(), DoubleEq(0));
}

// Test that one Second1D object is divided correctly by another.
TEST_F(ASecond1DClass, CorrectForDivision_X_Z) {
  double x = 1.2568;
  double z = -4.2157;
  Second1D X(x, 1, 0);
  Second1D Z(z, 0, 0);

  Second1D r = X / Z;

  ASSERT_THAT(r.value(), DoubleEq(x / z));
  ASSERT_THAT(r.first(), DoubleEq(1.0 / z));
  ASSERT_THAT(r.second(), DoubleEq(0));
}

// Tests that the function 2*x^2 is correctly divided by the function 3*x^3.
TEST_F(ASecond1DClass, CorrectForDivision_2X2_3X3) {
  double x = 1.2568;
  Second1D X2(2 * x * x, 4 * x, 4);
  Second1D X(3 * x * x * x, 9 * x * x, 18 * x);

  Second1D r = X2 / X;

  ASSERT_THAT(r.value(), DoubleEq(2.0 / (3 * x)));
  ASSERT_THAT(r.first(), DoubleNear(-2.0 / (3 * x * x), 1e-10));
  ASSERT_THAT(r.second(), DoubleNear(4.0 / (3 * x * x * x), 1e-10));
}

// Tests that the square root function is correctly implemented.
TEST_F(ASecond1DClass, CorrectForSqrt_2x3) {
  double x = 1.2568;
  Second1D X(x, 1, 0);

  Second1D R = sqrt(X * X * X * 2);

  double r = std::sqrt(x * x * x * 2);

  ASSERT_THAT(R.value(), DoubleEq(r));
  ASSERT_THAT(R.first(), DoubleEq(0.5 * 1 / r * 2 * x * x * 3));
  ASSERT_THAT(R.second(), DoubleEq(3 * x / (2 * r)));
}

// Tests that the cosine function is correctly implemented.
TEST_F(ASecond1DClass, CorrectForCos_2x_sqrt_x2) {
  double x = 1.5682;
  Second1D X(x, 1, 0);

  Second1D R = cos(2 * X * sqrt(X * X));

  ASSERT_THAT(R.value(), DoubleEq(std::cos(2 * x * x)));
  ASSERT_THAT(R.first(), DoubleEq(-4 * x * std::sin(2 * x * x)));
  ASSERT_THAT(R.second(), DoubleEq(-4 * (std::cos(2 * x * x) * 4 * x * x + std::sin(2 * x * x))));
}

// Tests that the exponential function is correctly implemented.
TEST_F(ASecond1DClass, CorrectForExp_x2) {
  double x = 1.2568;
  Second1D X(x, 1, 0);

  Second1D E = exp(X * X);
  double expV = std::exp(x * x);

  ASSERT_THAT(E.value(), DoubleEq(expV));
  ASSERT_THAT(E.first(), DoubleEq(expV * 2 * x));
  ASSERT_THAT(E.second(), DoubleEq(expV * (2 * x * 2 * x + 2)));
}

// Tests that a division f/x^2 gives the correct result.
TEST_F(ASecond1DClass, CorrectForDivision) {
  double f = -1.111;
  double x = 1.2568;
  double x2 = x * x;

  Second1D X(x, 1, 0);
  auto result = f / (X * X);

  double expectedValue = f / x2;
  double expectedFirstDerivative = -2 * f / (x * x2);
  double expectedSecondDerivative = 6 * f / (x2 * x2);

  ASSERT_THAT(result.value(), DoubleEq(expectedValue));
  ASSERT_THAT(result.first(), DoubleEq(expectedFirstDerivative));
  ASSERT_THAT(result.second(), DoubleEq(expectedSecondDerivative));
}

// Tests that the division f/x^4 gives the correct result.
TEST_F(ASecond1DClass, CorrectForOtherDivision) {
  double f = -1.111;
  double x = 1.2568;
  double x2 = x * x;
  double x4 = x2 * x2;

  Second1D X(x, 1, 0);
  auto result = f / (X * X * X * X);

  double expectedValue = f / x4;
  double expectedFirstDerivative = -4 * f / (x4 * x);
  double expectedSecondDerivative = 20 * f / (x4 * x2);

  ASSERT_THAT(result.value(), DoubleEq(expectedValue));
  ASSERT_THAT(result.first(), DoubleEq(expectedFirstDerivative));
  ASSERT_THAT(result.second(), DoubleEq(expectedSecondDerivative));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
