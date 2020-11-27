/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/ExternalQC/Orca/OrcaHessianOutputParser.h>
#include <Utils/GeometricDerivatives/HessianUtilities.h>
#include <Utils/Geometry/ElementTypes.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <gmock/gmock.h>
#include <algorithm>
#include <boost/dll/runtime_symbol_info.hpp>

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
    auto pathToResources = boost::dll::program_location().parent_path();
    pathToResources /= "Resources";
    ExternalQC::OrcaHessianOutputParser arbitraryHessianParser((pathToResources / "orca_test_calc.hess").string());
    arbitraryHessian = arbitraryHessianParser.getHessian();
  }
};

TEST_F(AHessianUtilitiesTest, CanComputeMassweightedEigenvalues) {
  HessianUtilities mwHessianUtilities(arbitraryHessian, arbitraryElements, arbitraryPositions, true);
  // Get eigenvalues
  Eigen::Vector3d eigenvalues;
  eigenvalues = mwHessianUtilities.getInternalEigenvalues();
  // Convert to cm-1
  double wavenumberConversion = Constants::invCentimeter_per_hartree * sqrt(Constants::u_per_electronRestMass);
  Eigen::Array3d wavenumbers = wavenumberConversion * eigenvalues.array().sqrt();
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

TEST_F(AHessianUtilitiesTest, ThrowsEmptyInternalHessianException) {
  std::stringstream stream("1\n\n"
                           "H     0.0000000000    0.0000000000    0.0000000000\n");
  AtomCollection monoatomicAtomCollection = XyzStreamHandler::read(stream);
  PositionCollection position = monoatomicAtomCollection.getPositions();
  ElementTypeCollection element = monoatomicAtomCollection.getElements();
  Eigen::Matrix<double, 3, 3> threedimensionalHessian;
  threedimensionalHessian << 1.0, 2.0, 3.0, 2.0, 2.0, 4.0, 3.0, 4.0, 3.0;
  HessianUtilities mwHessianUtilities(threedimensionalHessian, element, position, true);
  ASSERT_THROW(mwHessianUtilities.getInternalEigenvalues(), EmptyInternalHessianException);
}

} // namespace Utils
} // namespace Scine
