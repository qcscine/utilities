/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/MolecularTrajectory.h>
#include <gmock/gmock.h>
#include <cmath>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

class ANormalModesTest : public Test {
 public:
  DisplacementCollection displ;
  PositionCollection pos;
  ElementTypeCollection elements;

 protected:
  void SetUp() override {
    elements.push_back(ElementType::none);
    elements.push_back(ElementType::none);
    elements.push_back(ElementType::none);
    elements.push_back(ElementType::none);

    displ.resize(4, 3);
    pos.resize(4, 3);

    pos.row(0) = Position(-6.0, 2.0, 7.0);
    pos.row(1) = Position(-6.0, 1.0, 6.0);
    pos.row(2) = Position(-6.0, 0.0, 5.0);
    pos.row(3) = Position(-6.0, -1.0, 4.0);

    displ.row(0) = Position(0.5, 0.0, 0.0);
    displ.row(1) = Position(0.0, 0.32, 0.28);
    displ.row(2) = Position(0.42, 0.0, 0.0);
    displ.row(3) = Position(0.1, 0.0, 0.7);
  }
};

// Tests the class NormalModesContainer
TEST_F(ANormalModesTest, WrongIndicesAreIdentified) {
  AtomCollection structure(elements, pos);

  double wavenumber = 42.42;
  NormalMode normalMode(wavenumber, displ);
  NormalModesContainer modesContainer;
  modesContainer.add(normalMode);

  Eigen::MatrixXd normalModes = modesContainer.getNormalModes();
  Eigen::VectorXd shouldBeMode = Eigen::Map<const Eigen::VectorXd>(displ.data(), displ.size());
  for (int i = 0; i < displ.size(); ++i) {
    ASSERT_DOUBLE_EQ(normalModes.col(0)(i), shouldBeMode(i));
  }

  EXPECT_THROW(modesContainer.getMode(-2), std::runtime_error);
  EXPECT_THROW(modesContainer.getMode(10), std::runtime_error);
}

TEST_F(ANormalModesTest, ReturnsCorrectModeAndMolecularTrajectoryOfMode) {
  AtomCollection structure(elements, pos);

  double wavenumber = 42.42;
  NormalMode normalMode(wavenumber, displ);
  NormalModesContainer modesContainer;
  modesContainer.add(normalMode);

  auto wavenumbers = modesContainer.getWaveNumbers();
  EXPECT_THAT(wavenumber, Eq(wavenumbers[0]));

  auto returnedMode = modesContainer.getMode(0);
  EXPECT_THAT(returnedMode, Eq(displ));

  double scalingFactor = 2.5;
  auto mt = modesContainer.getModeAsMolecularTrajectory(0, structure, scalingFactor);

  auto firstPositions = mt[0];
  auto centerPositions = mt[mt.size() / 2];
  auto lastPositions = mt[mt.size() - 1];
  auto lastPositionsRef = pos + displ * scalingFactor * sin(2 * M_PI * (mt.size() - 1) / mt.size());

  for (int i = 0; i < firstPositions.rows(); ++i) {
    for (int j = 0; j < firstPositions.cols(); ++j) {
      ASSERT_THAT(firstPositions(i, j), DoubleNear(pos(i, j), 1e-8));
      ASSERT_THAT(centerPositions(i, j), DoubleNear(pos(i, j), 1e-8));
      ASSERT_THAT(lastPositions(i, j), DoubleNear(lastPositionsRef(i, j), 1e-8));
    }
  }
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
