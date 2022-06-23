/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/GeometricDerivatives/NormalModeAnalysis.h>
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

// Tests consistency between calculation of MW/non-MW normal-modes
// The normal modes matrix was calculated with PySCF with HF/STO-3G
TEST_F(ANormalModesTest, CheckCorrectNormalizationOfNormalModes) {
  // Elements type
  ElementTypeCollection elements;
  elements.push_back(ElementType::H);
  elements.push_back(ElementType::O);
  elements.push_back(ElementType::H);

  // Positions
  PositionCollection positions;
  positions.resize(3, 3);
  positions.row(0) = Position(0.7493682, 0.0, 0.4424329);
  positions.row(1) = Position(0.0000000, 0.0, -0.1653507);
  positions.row(2) = Position(-0.7493682, 0.0, 0.4424329);

  // Hessian matrix
  Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(9, 9);
  hessian(0, 0) = 0.438912E+00;
  hessian(0, 2) = 0.276724E+00;
  hessian(2, 0) = 0.276724E+00;
  hessian(0, 3) = -0.401975E+00;
  hessian(3, 0) = -0.401975E+00;
  hessian(0, 5) = -0.216312E+00;
  hessian(5, 0) = -0.216312E+00;
  hessian(0, 6) = -0.369362E-01;
  hessian(6, 0) = -0.369362E-01;
  hessian(0, 8) = -0.604124E-01;
  hessian(8, 0) = -0.604124E-01;
  hessian(2, 2) = 0.300103E+00;
  hessian(2, 3) = -0.337137E+00;
  hessian(3, 2) = -0.337137E+00;
  hessian(2, 5) = -0.317450E+00;
  hessian(5, 2) = -0.317450E+00;
  hessian(2, 6) = 0.604124E-01;
  hessian(6, 2) = 0.604124E-01;
  hessian(2, 8) = 0.173467E-01;
  hessian(8, 2) = 0.173467E-01;
  hessian(3, 3) = 0.803951E+00;
  hessian(3, 6) = -0.401975E+00;
  hessian(6, 3) = -0.401975E+00;
  hessian(3, 8) = 0.337137E+00;
  hessian(8, 3) = 0.337137E+00;
  hessian(5, 5) = 0.634900E+00;
  hessian(5, 6) = 0.216312E+00;
  hessian(6, 5) = 0.216312E+00;
  hessian(5, 8) = -0.317450E+00;
  hessian(8, 5) = -0.317450E+00;
  hessian(6, 6) = 0.438912E+00;
  hessian(6, 8) = -0.276724E+00;
  hessian(8, 6) = -0.276724E+00;
  hessian(8, 8) = 0.300103E+00;

  // Tests the normalization
  auto normalizedContainer = NormalModeAnalysis::calculateNormalModes(hessian, elements, positions);
  auto unnormalizedContainer = NormalModeAnalysis::calculateNormalModes(hessian, elements, positions, false);
  for (int iMode = 0; iMode < 3; iMode++) {
    double modeNorm = Eigen::Map<Eigen::VectorXd>(const_cast<double*>(normalizedContainer.getMode(iMode).data()), 9).norm();
    ASSERT_THAT(modeNorm, DoubleNear(1., 1e-8));
  }

  // Check the coherence between the mw and non-mw normal-modes
  Eigen::VectorXd mwFromGaussian(3);
  mwFromGaussian << 1.0785, 1.0491, 1.0774;
  for (int iMode = 0; iMode < 3; iMode++) {
    double reducedMass =
        Eigen::Map<Eigen::VectorXd>(const_cast<double*>(unnormalizedContainer.getMode(iMode).data()), 9).norm();
    // We don't expect perfect equality since the management of the rotational modes is different
    ASSERT_THAT(1. / std::pow(reducedMass, 2), DoubleNear(mwFromGaussian[iMode], 2e-3));
  }
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
