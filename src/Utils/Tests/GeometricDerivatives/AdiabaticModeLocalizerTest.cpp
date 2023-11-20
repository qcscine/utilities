/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Constants.h"
#include "Utils/GeometricDerivatives/AdiabaticModeLocalizer.h"
#include "Utils/IO/ChemicalFileFormats/XyzStreamHandler.h"
#include "Utils/MSVCCompatibility.h"
#include "Utils/MolecularTrajectory.h"
#include <gmock/gmock.h>
#include <Eigen/Eigen>
#include <algorithm>

using namespace testing;

namespace Scine {
namespace Utils {

class AnAdiabaticModeLocalizerTest : public Test {
 protected:
  // Input data
  Utils::AtomCollection atoms;
  Utils::BondOrderCollection bondOrders;
  std::vector<std::pair<int, int>> bonds;
  Eigen::MatrixXd hessian;

  // Reference results calculated with Compliance 3.0.2
  // http://www.oc.tu-bs.de/Grunenberg/compliance.html
  // Brandhorst, K.; Grunenberg, J.; Chem. Soc. Rev. 2008, 37 (8), 1558–1567. https://doi.org/10.1039/B717781J.
  // Brandhorst, K.; Grunenberg, J.; J. Chem. Phys. 2010, 132 (18), 184101. https://doi.org/10.1063/1.3413528.
  // Diagonal elements of compliance matrix for stretching movement in cm/N
  double compliance01 = 0.129; // O-H
  double compliance02 = 0.129; // O-H
  double compliance12 = 0.618; // H-H
  double unitConversion = 0.01 * Constants::bohr_per_meter * Constants::bohr_per_meter *
                          Constants::joule_per_hartree; // conversion from cm/N to atomic units i.e. bohr^2 / Hartree
  Eigen::SparseMatrix<double> forceConstants;

  // Corresponding localized cartesian modes
  std::map<std::pair<int, int>, DisplacementCollection> modes;
  DisplacementCollection mode01, mode02, mode12;

 private:
  void SetUp() final {
    std::stringstream stream("3\n\n"
                             "O   0.00000000000000      0.00000000006948      0.60066914031326\n"
                             "H   0.00000000000000      0.75798277245186     -0.01183457012843\n"
                             "H   0.00000000000000     -0.75798277252134     -0.0118345701848\n");

    atoms = XyzStreamHandler::read(stream);
    bondOrders.resize(3);
    bondOrders.setOrder(0, 1, 1.0);        // O-H
    bondOrders.setOrder(0, 2, 1.0);        // O-H
    bondOrders.setOrder(1, 2, 1.0);        // H-H
    bonds.push_back(std::make_pair(1, 2)); // H-H
    hessian.resize(9, 9);
    hessian << -4.56870e-04, 0.00000e+00, 0.00000e+00, 2.28430e-04, 0.00000e+00, 0.00000e+00, 2.28430e-04, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 6.16840e-01, 0.00000e+00, 0.00000e+00, -3.08420e-01, 2.49450e-01, 0.00000e+00,
        -3.08420e-01, -2.49450e-01, 0.00000e+00, 0.00000e+00, 4.38716e-01, 0.00000e+00, 1.89848e-01, -2.19358e-01,
        0.00000e+00, -1.89848e-01, -2.19358e-01, 2.28430e-04, 0.00000e+00, 0.00000e+00, -1.75460e-04, 0.00000e+00,
        0.00000e+00, -5.29700e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.08420e-01, 1.89848e-01, 0.00000e+00,
        3.40239e-01, -2.19649e-01, 0.00000e+00, -3.18191e-02, 2.98012e-02, 0.00000e+00, 2.49450e-01, -2.19358e-01,
        0.00000e+00, -2.19649e-01, 2.10526e-01, 0.00000e+00, -2.98012e-02, 8.83245e-03, 2.28430e-04, 0.00000e+00,
        0.00000e+00, -5.29700e-05, 0.00000e+00, 0.00000e+00, -1.75460e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        -3.08420e-01, -1.89848e-01, 0.00000e+00, -3.18191e-02, -2.98012e-02, 0.00000e+00, 3.40239e-01, 2.19649e-01,
        0.00000e+00, -2.49450e-01, -2.19358e-01, 0.00000e+00, 2.98012e-02, 8.83245e-03, 0.00000e+00, 2.19649e-01,
        2.10526e-01;

    // Convert compliance matrix elements in cm/N to relaxed force constants in atomic units
    forceConstants.resize(3, 3);
    forceConstants.coeffRef(0, 1) = 1 / (compliance01 * unitConversion);
    forceConstants.coeffRef(0, 2) = 1 / (compliance02 * unitConversion);
    forceConstants.coeffRef(1, 2) = 1 / (compliance12 * unitConversion);
    forceConstants = forceConstants + Eigen::SparseMatrix<double>(forceConstants.transpose());
    // Reference displacements calculated with Compliance 3.0.2
    mode01.resize(3, 3);
    mode01 << -7.69581e-17, -0.690081, 0.630369, 2.68039e-17, 0.907821, -0.594148, 5.87917e-17, -0.21774, -0.036221;
    mode02.resize(3, 3);
    mode02 << 8.64258e-17, 0.690081, 0.630369, 9.46014e-17, 0.21774, -0.036221, -1.66067e-16, -0.907821, -0.594148;
    mode12.resize(3, 3);
    mode12 << 2.02193e-16, 1.07692e-14, -2.77727, -1.57597e-15, 4.81346, 1.38864, 2.10847e-15, -4.81346, 1.38864;
    modes.insert(std::pair<std::pair<int, int>, DisplacementCollection>(std::make_pair(0, 1), mode01));
    modes.insert(std::pair<std::pair<int, int>, DisplacementCollection>(std::make_pair(0, 2), mode02));
    modes.insert(std::pair<std::pair<int, int>, DisplacementCollection>(std::make_pair(1, 2), mode12));
  }
};

TEST_F(AnAdiabaticModeLocalizerTest, CanStoreInAdiabaticModesContainer) {
  // Test the basic functionality of the modes container
  AdiabaticModesContainer modesContainer;
  double wavenumber = 42.42;
  double forceConstant = 24.42;
  NormalMode normalMode(wavenumber, mode01);
  std::pair<int, int> bond(2, 0);
  std::pair<int, int> reverseBond(0, 2);
  modesContainer.addMode(bond, normalMode, forceConstant);

  // Throws exception if trying to overwrite
  EXPECT_THROW(modesContainer.addMode(bond, normalMode, forceConstant), AdiabaticLocalizationExceptions::modeReadded);
  EXPECT_THROW(modesContainer.addMode(reverseBond, normalMode, forceConstant), AdiabaticLocalizationExceptions::modeReadded);

  // Gets wavenumber
  auto readWavenumbers = modesContainer.getWaveNumbers();
  EXPECT_THAT(readWavenumbers.find(bond)->second, Eq(wavenumber));
  EXPECT_THAT(readWavenumbers.find(reverseBond)->second, Eq(wavenumber));

  // Gets force constants
  auto readForceConstants = modesContainer.getForceConstants();
  EXPECT_THAT(readForceConstants.find(bond)->second, Eq(forceConstant));
  EXPECT_THAT(readForceConstants.find(reverseBond)->second, Eq(forceConstant));

  // Gets mode
  auto returnedMode = modesContainer.getMode(bond);
  EXPECT_THAT(returnedMode, Eq(mode01));
  returnedMode = modesContainer.getMode(reverseBond);
  EXPECT_THAT(returnedMode, Eq(mode01));

  // Gets molecular trajectory
  double scalingFactor = 2.5;
  MolecularTrajectory trajectory = modesContainer.getModeAsMolecularTrajectory(reverseBond, atoms, scalingFactor);

  auto positions = atoms.getPositions();
  auto firstPositions = trajectory[0];
  auto centerPositions = trajectory[trajectory.size() / 2];
  auto lastPositions = trajectory[trajectory.size() - 1];
  auto lastPositionsRef = positions + mode01 * scalingFactor * sin(2 * M_PI * (trajectory.size() - 1) / trajectory.size());

  for (int i = 0; i < firstPositions.rows(); ++i) {
    for (int j = 0; j < firstPositions.cols(); ++j) {
      ASSERT_THAT(firstPositions(i, j), DoubleNear(positions(i, j), 1e-8));
      ASSERT_THAT(centerPositions(i, j), DoubleNear(positions(i, j), 1e-8));
      ASSERT_THAT(lastPositions(i, j), DoubleNear(lastPositionsRef(i, j), 1e-8));
    }
  }
}

TEST_F(AnAdiabaticModeLocalizerTest, CanConstructModeLocalizerFromBondOrders) {
  AdiabaticModeLocalizer modeLocalizer(hessian, atoms, bondOrders);
}

TEST_F(AnAdiabaticModeLocalizerTest, CanConstructModeLocalizerFromBonds) {
  AdiabaticModeLocalizer modeLocalizer(hessian, atoms, bonds);
}

TEST_F(AnAdiabaticModeLocalizerTest, CorrectlyLocalizesWaterModes) {
  AdiabaticModeLocalizer modeLocalizer(hessian, atoms, bondOrders);
  AdiabaticModesContainer calculatedAdiabaticModes = modeLocalizer.localizeModes();
  // Check force constants
  std::map<std::pair<int, int>, double> calculatedForceConstants = calculatedAdiabaticModes.getForceConstants();
  for (auto const& it : calculatedForceConstants) {
    ASSERT_THAT(it.second, DoubleNear(forceConstants.coeff(it.first.first, it.first.second), 2e-3));
  }
  // Check modes
  for (auto const& it : calculatedForceConstants) {
    // Get normalized modes
    DisplacementCollection calculatedMode = calculatedAdiabaticModes.getMode(it.first);
    calculatedMode.normalize();
    DisplacementCollection referenceMode = modes.find(std::minmax(it.first.first, it.first.second))->second;
    referenceMode.normalize();

    ASSERT_TRUE(calculatedMode.isApprox(referenceMode, 1e-2));
  }
}

// For a diatomic molecule the adiabatic mode localizer should reproduce the vibrational normal mode
TEST_F(AnAdiabaticModeLocalizerTest, CanReproduceDiatomicNormalMode) {
  // O2
  std::stringstream stream("2\n\n"
                           "O   0.00000000000000      0.00000000000000      1.10808053007267\n"
                           "O   0.00000000000000      0.00000000000000     -0.10808053007267\n");
  AtomCollection o2Atoms = XyzStreamHandler::read(stream);
  // Hessian calculated with Orca/4.2.0
  Eigen::MatrixXd o2Hessian(6, 6);
  o2Hessian << -7.4989979740E-05, 4.5052586967E-07, 6.4127257261E-11, 7.4989979113E-05, -4.5052596080E-07, -6.4028648945E-11,
      4.5235939859E-07, -2.1195184661E-04, 6.1101522459E-11, -4.5235913895E-07, 2.1195184783E-04, -6.1093416715E-11,
      6.2822856365E-11, 5.8773212750E-11, 7.6393436568E-01, -6.3587351999E-11, -5.8916513479E-11, -7.6393436568E-01,
      7.4989978485E-05, -4.5051780426E-07, 1.0507663671E-11, -7.4989979312E-05, 4.5051804354E-07, -1.1363319352E-11,
      -4.5235045930E-07, 2.1195182889E-04, -6.4914792451E-12, 4.5235054399E-07, -2.1195182763E-04, 6.4308882910E-12,
      8.4502591903E-12, -8.1052935718E-12, -7.6393436583E-01, -8.4320950148E-12, 8.2175919268E-12, 7.6393436583E-01;

  BondOrderCollection o2BondOrders(2);
  o2BondOrders.setOrder(0, 1, 1.0);

  // Reference data calculated with Orca/4.2.0
  double wavenumber = 1588.52; // in cm-1
  DisplacementCollection normalMode;
  normalMode.resize(2, 3);
  // Normalized mode
  normalMode << 0.000000, 0.000000, 0.707107, 0.000000, 0.000000, -0.707107;

  // Localize to bond
  AdiabaticModeLocalizer modeLocalizer(o2Hessian, o2Atoms, o2BondOrders);
  AdiabaticModesContainer adiabaticModes = modeLocalizer.localizeModes();
  // Check localized mode
  DisplacementCollection adiabaticMode = adiabaticModes.getMode(std::make_pair(0, 1));
  adiabaticMode.normalize();
  ASSERT_TRUE(adiabaticMode.isApprox(normalMode, 1e-2));
  // Check wavenumber
  std::map<std::pair<int, int>, double> adiabaticWavenumbers = adiabaticModes.getWaveNumbers();
  ASSERT_THAT(adiabaticWavenumbers.find(std::make_pair(0, 1))->second, DoubleNear(wavenumber, 1e-1));
}

} // namespace Utils
} // namespace Scine
