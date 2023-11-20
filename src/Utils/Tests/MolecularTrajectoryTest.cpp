/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/MolecularTrajectory.h"
#include "Utils/Typenames.h"
#include "gmock/gmock.h"

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class MolecularTrajectoryTest MolecularTrajectoryTest.cpp
 * @brief Comprises separate tests for the MolecularTrajectory class,
 * which is mostly tested in the MolecularDynamicsTest
 * @test
 */
class MolecularTrajectoryTest : public Test {
 public:
  ElementTypeCollection elements;
  PositionCollection p1, p2, p3;

 protected:
  void SetUp() override {
    elements = {ElementType::H, ElementType::F};
    p1 = Eigen::MatrixX3d::Zero(2, 3);
    p2 = Eigen::MatrixX3d::Zero(2, 3);
    p3 = Eigen::MatrixX3d::Zero(2, 3);
    // clang-format off
    p1 <<  0.0, 1.0, 2.0,
           0.5, 1.5, 2.5;
    p2 <<  1.0, 2.0, 3.0,
           1.5, 2.5, 3.5;
    p3 << 10.0, 2.0, 3.0,
           1.5, 2.5, 3.5;
    // clang-format on
  }
};

TEST_F(MolecularTrajectoryTest, RmsdLimitWorks) {
  auto traj = MolecularTrajectory(elements, 1.5);
  ASSERT_TRUE(traj.getElementTypes().size() == 2);
  ASSERT_TRUE(traj.getElementTypes()[0] == ElementType::H);
  ASSERT_TRUE(traj.getElementTypes()[1] == ElementType::F);
  traj.push_back(p1);
  traj.push_back(p2);
  traj.push_back(p2);
  traj.push_back(p2);
  traj.push_back(p1);
  traj.push_back(p3);
  ASSERT_TRUE(traj.size() == 4);
  auto stricter_traj = MolecularTrajectory(5.0);
  stricter_traj.setElementTypes(elements);
  stricter_traj.push_back(p1);
  stricter_traj.push_back(p2);
  stricter_traj.push_back(p2);
  stricter_traj.push_back(p2);
  stricter_traj.push_back(p1);
  ASSERT_TRUE(stricter_traj.size() == 1);
  stricter_traj.push_back(p3);
  ASSERT_TRUE(stricter_traj.size() == 2);
  stricter_traj.push_back(p1);
  ASSERT_TRUE(stricter_traj.size() == 3);
  stricter_traj.push_back(p1);
  ASSERT_TRUE(stricter_traj.size() == 3);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
