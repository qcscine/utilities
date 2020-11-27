/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Utils/Math/AutomaticDifferentiation/Second1D.h>
#include <Utils/Math/FullSecondDerivativeCollection.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
using namespace AutomaticDifferentiation;
namespace Tests {

// Test class DerivativesContainerTest
// This class tests the FullSecondDerivativeCollection class.
class DerivativesContainerTest : public Test {
 public:
 protected:
  void SetUp() override {
  }
};

// Tests that hessian and gradients are zero after call of setZero function.
TEST_F(DerivativesContainerTest, HessianAndGradientsAreZero) {
  int N = 5;
  FullSecondDerivativeCollection collection(N);
  collection.setZero();

  // get hessian and gradients from collection
  auto hessian = collection.getHessianMatrix();
  auto gradients = collection.getReferenceGradients();

  // assert that all elements are zero
  for (int i = 0; i < 3 * N; ++i) {
    for (int j = 0; j < 3; ++j) {
      ASSERT_THAT(gradients(i / 3, j), 0);
    }
    for (int k = 0; k < N; ++k) {
      ASSERT_THAT(hessian(i, k), 0);
    }
  }
}

// Tests that a derivative is correctly added to the container.
TEST_F(DerivativesContainerTest, AddDerivativeCorrectlyToContainer) {
  int N = 2;
  FullSecondDerivativeCollection collection(N);
  collection.setZero();

  // Define diatomic molecule. Define vector r between the two atoms as Second1D object
  Position posA(0.0, 0.0, 0.0);
  Position posB(0.0, 0.0, 1.5);
  Eigen::Vector3d rVector = posB - posA;
  Second1D r(rVector.norm(), 1, 0);

  // Assume, E = 5*r^2
  auto energy = 5 * r * r;
  Second3D v = get3Dfrom1D<DerivativeOrder::Two>(energy, rVector);

  // Function to be tested
  collection.addDerivative(0, 1, v);

  // get hessian and gradients from collection
  auto hessian = collection.getHessianMatrix();
  auto gradients = collection.getReferenceGradients();

  // Check that the gradients are correct
  ASSERT_THAT(gradients(0, 0), 0);
  ASSERT_THAT(gradients(0, 1), 0);
  ASSERT_THAT(gradients(0, 2), -15);
  ASSERT_THAT(gradients(1, 0), 0);
  ASSERT_THAT(gradients(1, 1), 0);
  ASSERT_THAT(gradients(1, 2), 15);

  // Check some hessian elements for correctness, too
  ASSERT_THAT(hessian(0, 0), 10);
  ASSERT_THAT(hessian(0, 1), 0);
  ASSERT_THAT(hessian(0, 2), 0);
  ASSERT_THAT(hessian(0, 3), -10);

  // Check that hessian is symmetrical
  for (int i = 1; i < 3 * N; ++i) {
    for (int j = 0; j < i; ++j) {
      ASSERT_THAT(hessian(i, j), hessian(j, i));
    }
  }
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
