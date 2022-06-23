/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Math/LinearSumAssignment.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Utils {

using namespace LinearSumAssignment;
using namespace testing;

class ALinearSumAssignmentTest : public Test {};

TEST_F(ALinearSumAssignmentTest, CanAssignWithQuadraticIdentityScore) {
  Eigen::MatrixXd cost = Eigen::MatrixXd::Identity(20, 20);
  Solver<LinearSumAssignment::Type::Maximize> solver(cost);
  for (const auto it : solver()) {
    EXPECT_EQ(it.first, it.second);
  }
}

TEST_F(ALinearSumAssignmentTest, CanAssignWithQuadraticIdentityCost) {
  Eigen::MatrixXd cost = Eigen::MatrixXd::Ones(20, 20) - Eigen::MatrixXd::Identity(20, 20);
  Solver<LinearSumAssignment::Type::Minimize> solver(cost);
  for (const auto it : solver()) {
    EXPECT_EQ(it.first, it.second);
  }
}
TEST_F(ALinearSumAssignmentTest, CanAssignWithRectangularIdentityCost) {
  Solver<LinearSumAssignment::Type::Maximize> solver(Eigen::MatrixXd::Identity(20, 30));
  for (const auto it : solver()) {
    EXPECT_EQ(it.first, it.second);
  }
}

TEST_F(ALinearSumAssignmentTest, CanAssignWithRectangularIdentityCostTransposed) {
  Solver<LinearSumAssignment::Type::Maximize> solver(Eigen::MatrixXd::Identity(30, 20));
  for (const auto it : solver()) {
    EXPECT_EQ(it.first, it.second);
  }
}

TEST_F(ALinearSumAssignmentTest, CanAssignAsScipyImplementationTransposed) {
  Eigen::MatrixXd arbitraryCost(8, 8);

  arbitraryCost << 0.51154057, 1.74065748, 1.84011773, 1.65797869, 1.87969028, 1.04066876, 1.1764934, 1.43272137,
      1.52843714, 1.04591385, 0.86802574, 1.14905381, 0.72940933, 0.75420008, 1.7549163, 1.49694873, 0.88240342,
      1.41762155, 1.013285, 1.62462635, 1.05332804, 0.15128548, 1.92434222, 0.83774796, 1.79080658, 1.85844126, 0.52789256,
      1.16508171, 1.83252667, 0.42761326, 1.58522088, 0.77804918, 0.09692329, 1.2155551, 0.85463422, 1.77803774,
      0.88354351, 1.8790176, 1.18369512, 1.18370504, 1.63625258, 1.64351976, 1.10080942, 1.70957339, 1.58388997,
      1.41128758, 0.84302794, 1.84760905, 0.691641, 0.09446226, 1.97218777, 1.38822919, 1.50325824, 1.41653303, 0.15864116,
      0.14030929, 1.96456312, 0.58830223, 0.69833206, 0.48575604, 0.62108857, 1.80657139, 1.78996699, 0.5179649;

  arbitraryCost.transposeInPlace();
  Solver<Type::Minimize> solver(arbitraryCost);

  // From scipy.optimize.linear_sum_assignment()
  Eigen::VectorXi expected(8);
  expected << 4, 6, 3, 7, 1, 2, 5, 0;
  auto actual = solver();

  for (int i = 0; i < 8; ++i) {
    EXPECT_EQ(actual.at(i), expected(i));
  }

  EXPECT_NEAR(solver.cost(), 4.36147827, 1e-7);
}

TEST_F(ALinearSumAssignmentTest, CanAssignAsScipyImplementation) {
  Eigen::MatrixXd arbitraryCost(8, 8);

  arbitraryCost << 0.51154057, 1.74065748, 1.84011773, 1.65797869, 1.87969028, 1.04066876, 1.1764934, 1.43272137,
      1.52843714, 1.04591385, 0.86802574, 1.14905381, 0.72940933, 0.75420008, 1.7549163, 1.49694873, 0.88240342,
      1.41762155, 1.013285, 1.62462635, 1.05332804, 0.15128548, 1.92434222, 0.83774796, 1.79080658, 1.85844126, 0.52789256,
      1.16508171, 1.83252667, 0.42761326, 1.58522088, 0.77804918, 0.09692329, 1.2155551, 0.85463422, 1.77803774,
      0.88354351, 1.8790176, 1.18369512, 1.18370504, 1.63625258, 1.64351976, 1.10080942, 1.70957339, 1.58388997,
      1.41128758, 0.84302794, 1.84760905, 0.691641, 0.09446226, 1.97218777, 1.38822919, 1.50325824, 1.41653303, 0.15864116,
      0.14030929, 1.96456312, 0.58830223, 0.69833206, 0.48575604, 0.62108857, 1.80657139, 1.78996699, 0.5179649;

  Solver<Type::Minimize> solver(arbitraryCost);

  // From scipy.optimize.linear_sum_assignment()
  Eigen::VectorXi expected(8);
  expected << 7, 4, 5, 2, 0, 6, 1, 3;
  auto actual = solver();

  for (int i = 0; i < 8; ++i) {
    EXPECT_EQ(actual.at(i), expected(i));
  }

  EXPECT_NEAR(solver.cost(), 4.36147827, 1e-7);
}
} // namespace Utils
} // namespace Scine
