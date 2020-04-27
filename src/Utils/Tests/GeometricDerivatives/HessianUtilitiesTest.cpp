/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Tests/ExternalQC/externalQC_output_location.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/ExternalQC/Orca/OrcaHessianOutputParser.h>
#include <Utils/GeometricDerivatives/HessianUtilities.h>
#include <Utils/Geometry/ElementTypes.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <gmock/gmock.h>
#include <algorithm>

using namespace testing;
namespace Scine {
namespace Utils {

class AHessianUtilitiesTest : public Test {
 public:
  ElementTypeCollection arbitraryElements;
  PositionCollection arbitraryPositions;
  HessianMatrix arbitraryHessian;

 private:
  void SetUp() final {
    // Use the ORCA output available in the ExternalQC tests
    // Generate AtomCollection corresponding to H2O+ as given in ORCA input
    std::stringstream stream("3\n\n"
                             "H      0.7493682000    0.0000000000    0.4424329000\n"
                             "O      0.0000000000    0.0000000000   -0.1653507000\n"
                             "H     -0.7493682000    0.0000000000    0.4424329000\n");
    AtomCollection arbitraryAtomCollection = XyzStreamHandler::read(stream);
    arbitraryPositions = arbitraryAtomCollection.getPositions();
    arbitraryElements = arbitraryAtomCollection.getElements();
    ExternalQC::OrcaHessianOutputParser arbitraryHessianParser(orca_test_hessian_output);
    arbitraryHessian = arbitraryHessianParser.getHessian();
  }
};

TEST_F(AHessianUtilitiesTest, CanComputeMassweightedEigenvalues) {
  HessianUtilities mwHessianUtilities(arbitraryHessian, arbitraryElements, arbitraryPositions, true);
  // Get eigenvalues
  Eigen::Vector3d eigenvalues;
  eigenvalues = mwHessianUtilities.getInternalEigenvalues();
  // Compute frequencies and transform to wavenumber
  double freqToWvn = 5140.48686; // Conversion factor between frequency in a.u. and wavenumber in cm^-1
  Eigen::Array3d wavenumbers = freqToWvn * eigenvalues.array().sqrt();
  // Compare to ORCA result; Slight deviations expected bc of different methods to project out rotation
  // EigenSolver returns sorted eigenvalues
  ASSERT_THAT(wavenumbers(0), DoubleNear(1348.49, 1));
  ASSERT_THAT(wavenumbers(1), DoubleNear(3929.92, 1));
  ASSERT_THAT(wavenumbers(2), DoubleNear(3944.83, 1));
}

TEST_F(AHessianUtilitiesTest, CanComputeMassweightedCartesianDisplacements) {
  HessianUtilities mwHessianUtilities(arbitraryHessian, arbitraryElements, arbitraryPositions, true);
  Eigen::Matrix<double, 9, 3> cartesianDisplacements = mwHessianUtilities.getBackTransformedInternalEigenvectors();
  // Orca normal modes
  Eigen::Matrix<double, 9, 3> orcaModes;
  // clang-format off
  orcaModes << -0.442359,   0.574440,   0.547872,
                0.000000,   0.000000,   0.000000,
                0.549474,   0.410706,   0.444358,
                0.000000,  -0.000000,  -0.069036,
                0.000000,   0.000000,   0.000000,
               -0.069238,  -0.051752,   0.000000,
                0.442359,  -0.574440,   0.547872,
                0.000000,   0.000000,   0.000000,
                0.549474,   0.410706,  -0.444358;
  // clang-format on

  // Compare to Orca normal modes
  // Eigenvectors are sorted according to corresponding eigenvalues
  // Modes are normalized but could still have the opposite sign
  for (unsigned int i = 0; i < 3; ++i) {
    ASSERT_THAT(std::min((orcaModes.col(i) - cartesianDisplacements.col(i)).squaredNorm(),
                         (orcaModes.col(i) + cartesianDisplacements.col(i)).squaredNorm()),
                DoubleNear(0, 1e-6));
  }
}

} // namespace Utils
} // namespace Scine
