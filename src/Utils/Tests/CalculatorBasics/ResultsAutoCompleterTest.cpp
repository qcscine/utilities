/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/ResultsAutoCompleter.h>
#include <Utils/ExternalQC/Orca/OrcaHessianOutputParser.h>
#include <Utils/Geometry/ElementTypes.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <gmock/gmock.h>
#include <Eigen/Eigenvalues>
#include <boost/dll/runtime_symbol_info.hpp>

using namespace testing;
namespace Scine {
namespace Utils {

class AResultsAutoCompleterTest : public Test {
 public:
  Results arbitraryResults;
  std::unique_ptr<ResultsAutoCompleter> arbitraryResultsAutoCompleter;
  HessianMatrix arbitraryHessian;
  Eigen::MatrixXd arbitraryOneElectronMatrix;
  BondOrderCollection arbitraryBondOrders;
  ElementTypeCollection arbitraryElements;
  PositionCollection randomPositions;
  AtomCollection randomAtomCollection;
  MolecularOrbitals arbitraryCoefficientMatrix;
  LcaoUtils::ElectronicOccupation arbitraryOccupation;
  Eigen::MatrixXd arbitraryOverlapMatrix;
  DensityMatrix arbitraryDensityMatrix;
  AtomsOrbitalsIndexes arbitraryAOToAtomMapping;
  Eigen::MatrixXd generateRandomOverlapMatrix(int nOrbitals);
  Eigen::MatrixXd generateNonOrthogonalCoefficientMatrix(const Eigen::MatrixXd& overlap, int nOrbitals);

 private:
  void SetUp() final {
    arbitraryElements = ElementTypeCollection{ElementType::Am, ElementType::Ce, ElementType::Ca};
    randomPositions = Eigen::MatrixX3d::Random(3, 3);
    randomAtomCollection.resize(arbitraryElements.size());
    randomAtomCollection.setElements(arbitraryElements);
    randomAtomCollection.setPositions(randomPositions);
    const int nOrbitals = 15;
    arbitraryOccupation.fillLowestRestrictedOrbitalsWithElectrons(4);
    arbitraryOverlapMatrix = generateRandomOverlapMatrix(nOrbitals);
    // Generate coefficient matrix as Eigen type
    auto arbitraryCoefficients = generateNonOrthogonalCoefficientMatrix(arbitraryOverlapMatrix, nOrbitals);
    // Deduce coefficients as MolecularOrbitals type
    arbitraryCoefficientMatrix = MolecularOrbitals::createFromRestrictedCoefficients(arbitraryCoefficients);
    // Set up AO to Atom Mapping
    arbitraryAOToAtomMapping.clear();
    arbitraryAOToAtomMapping.setSize(3);
    arbitraryAOToAtomMapping.addAtom(7);
    arbitraryAOToAtomMapping.addAtom(3);
    arbitraryAOToAtomMapping.addAtom(5);
  }
};

// Adapted from HFWaveFunctionOverlapTest to get arbitrary overlap and MO coefficients
Eigen::MatrixXd AResultsAutoCompleterTest::generateRandomOverlapMatrix(int nOrbitals) {
  Eigen::MatrixXd random = Eigen::MatrixXd::Random(nOrbitals, nOrbitals);
  random /= 10;
  for (int i = 0; i < nOrbitals; ++i)
    random(i, i) = 1;
  return .5 * (random + random.transpose());
}
Eigen::MatrixXd AResultsAutoCompleterTest::generateNonOrthogonalCoefficientMatrix(const Eigen::MatrixXd& overlap, int nOrbitals) {
  Eigen::MatrixXd matrixToDiagonalize = 20 * Eigen::MatrixXd::Random(nOrbitals, nOrbitals);
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es;
  es.compute(matrixToDiagonalize, overlap);
  return es.eigenvectors();
}

TEST_F(AResultsAutoCompleterTest, CanConstructResultAutoCompleter) {
  arbitraryResultsAutoCompleter = std::make_unique<ResultsAutoCompleter>(randomAtomCollection);
}

TEST_F(AResultsAutoCompleterTest, CanDetermineDensityMatrixGenerability) {
  arbitraryResultsAutoCompleter = std::make_unique<ResultsAutoCompleter>(randomAtomCollection);
  arbitraryResults.set<Property::SuccessfulCalculation>(true);
  arbitraryResults.set<Property::ElectronicOccupation>(arbitraryOccupation);
  // Just one more random property
  arbitraryResults.set<Property::OneElectronMatrix>(arbitraryOneElectronMatrix);

  ASSERT_FALSE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::DensityMatrix));

  arbitraryResults.set<Property::CoefficientMatrix>(arbitraryCoefficientMatrix);

  ASSERT_TRUE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::DensityMatrix));
}

TEST_F(AResultsAutoCompleterTest, CanDetermineThermochemistryGenerability) {
  arbitraryResultsAutoCompleter = std::make_unique<ResultsAutoCompleter>(randomAtomCollection);
  arbitraryResults.set<Property::SuccessfulCalculation>(true);
  arbitraryResults.set<Property::ElectronicOccupation>(arbitraryOccupation);
  // Just one more random property
  arbitraryResults.set<Property::OneElectronMatrix>(arbitraryOneElectronMatrix);

  ASSERT_FALSE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::Thermochemistry));

  arbitraryResults.set<Property::Hessian>(arbitraryHessian);

  ASSERT_TRUE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::Thermochemistry));
}

TEST_F(AResultsAutoCompleterTest, CanDetermineAtomicChargesGenerability) {
  arbitraryResultsAutoCompleter = std::make_unique<ResultsAutoCompleter>(randomAtomCollection);
  arbitraryResults.set<Property::SuccessfulCalculation>(true);
  arbitraryResults.set<Property::AOtoAtomMapping>(arbitraryAOToAtomMapping);
  arbitraryResults.set<Property::DensityMatrix>(arbitraryDensityMatrix);

  ASSERT_FALSE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::AtomicCharges));

  arbitraryResults.set<Property::OverlapMatrix>(arbitraryOverlapMatrix);

  ASSERT_TRUE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::AtomicCharges));
}

TEST_F(AResultsAutoCompleterTest, CanDetermineBondOrderGenerability) {
  arbitraryResultsAutoCompleter = std::make_unique<ResultsAutoCompleter>(randomAtomCollection);
  arbitraryResults.set<Property::SuccessfulCalculation>(true);
  arbitraryResults.set<Property::AOtoAtomMapping>(arbitraryAOToAtomMapping);
  arbitraryResults.set<Property::DensityMatrix>(arbitraryDensityMatrix);

  ASSERT_FALSE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::BondOrderMatrix));

  arbitraryResults.set<Property::OverlapMatrix>(arbitraryOverlapMatrix);

  ASSERT_TRUE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::BondOrderMatrix));
}

TEST_F(AResultsAutoCompleterTest, CanDetermineEnergyUngenerability) {
  arbitraryResultsAutoCompleter = std::make_unique<ResultsAutoCompleter>(randomAtomCollection);
  arbitraryResults.set<Property::SuccessfulCalculation>(true);
  ASSERT_FALSE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::Energy));
}

TEST_F(AResultsAutoCompleterTest, CanComputeThermochemistry) {
  // Generate AtomCollection corresponding to H2O+ as given in ORCA input
  std::stringstream stream("3\n\n"
                           "H      0.7493682000    0.0000000000    0.4424329000\n"
                           "O      0.0000000000    0.0000000000   -0.1653507000\n"
                           "H     -0.7493682000    0.0000000000    0.4424329000\n");
  auto structure = XyzStreamHandler::read(stream);
  arbitraryResultsAutoCompleter = std::make_unique<ResultsAutoCompleter>(structure);
  // C2v symmetry
  arbitraryResultsAutoCompleter->setMolecularSymmetryNumber(2);

  arbitraryResults.set<Property::SuccessfulCalculation>(true);
  // Set electronic occupation
  // Only thing that is really needed from occupation is multiplicity so make sure it is correct
  LcaoUtils::ElectronicOccupation unrestrictedOccupation;
  unrestrictedOccupation.fillLowestUnrestrictedOrbitals(5, 4);
  arbitraryResults.set<Property::ElectronicOccupation>(unrestrictedOccupation);
  // Electronic energy from Resources/orca_test_calc.out
  arbitraryResults.set<Property::Energy>(-75.818269296087);

  // Read in Hessian from ORCA output given in ExternalQC test
  auto pathToResources = boost::dll::program_location().parent_path();
  pathToResources /= "Resources";
  const auto hessianFile = (pathToResources / "orca_test_calc.hess").string();
  HessianMatrix hessian = ExternalQC::OrcaHessianOutputParser::getHessian(hessianFile);
  arbitraryResults.set<Property::Hessian>(hessian);

  // Should be able to compute Thermochemistry
  ASSERT_TRUE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::Thermochemistry));
  // Launch calculation
  arbitraryResultsAutoCompleter->generateProperties(arbitraryResults, structure);

  // Should have computed thermochemical properties
  ASSERT_TRUE(arbitraryResults.has<Property::Thermochemistry>());
  // Insufficient information for calculating Atomic Charges so they should not be there
  ASSERT_FALSE(arbitraryResults.has<Property::AtomicCharges>());

  // Check whether computed Thermochemical properties match orca_test_calc.out
  ASSERT_THAT(arbitraryResults.get<Property::Thermochemistry>().vibrationalComponent.zeroPointVibrationalEnergy,
              DoubleNear(0.02101209, 1e-5));
  ASSERT_THAT(arbitraryResults.get<Property::Thermochemistry>().overall.enthalpy, DoubleNear(-75.79347127, 1e-5));
  ASSERT_THAT(arbitraryResults.get<Property::Thermochemistry>().overall.entropy * 298.15, DoubleNear(0.02210368, 1e-5));
  ASSERT_THAT(arbitraryResults.get<Property::Thermochemistry>().overall.gibbsFreeEnergy, DoubleNear(-75.81557495, 1e-5));
}

TEST_F(AResultsAutoCompleterTest, CanComputeThermochemistryReverseOccupation) {
  // Generate AtomCollection corresponding to H2O+ as given in ORCA input
  std::stringstream stream("3\n\n"
                           "H      0.7493682000    0.0000000000    0.4424329000\n"
                           "O      0.0000000000    0.0000000000   -0.1653507000\n"
                           "H     -0.7493682000    0.0000000000    0.4424329000\n");
  auto structure = XyzStreamHandler::read(stream);
  arbitraryResultsAutoCompleter = std::make_unique<ResultsAutoCompleter>(structure);
  // C2v symmetry
  arbitraryResultsAutoCompleter->setMolecularSymmetryNumber(2);

  arbitraryResults.set<Property::SuccessfulCalculation>(true);
  // Set electronic occupation
  // Only thing that is really needed from occupation is multiplicity so make sure it is correct
  LcaoUtils::ElectronicOccupation unrestrictedOccupation;
  // Reverse occupation w.r.t. CanComputeThermochemistry
  unrestrictedOccupation.fillLowestUnrestrictedOrbitals(4, 5);
  arbitraryResults.set<Property::ElectronicOccupation>(unrestrictedOccupation);
  // Electronic energy from Resources/orca_test_calc.out
  arbitraryResults.set<Property::Energy>(-75.818269296087);

  // Read in Hessian from ORCA output given in ExternalQC test
  auto pathToResources = boost::dll::program_location().parent_path();
  pathToResources /= "Resources";
  const auto hessianFile = (pathToResources / "orca_test_calc.hess").string();
  HessianMatrix hessian = ExternalQC::OrcaHessianOutputParser::getHessian(hessianFile);
  arbitraryResults.set<Property::Hessian>(hessian);

  // Should be able to compute Thermochemistry
  ASSERT_TRUE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::Thermochemistry));
  // Launch calculation
  arbitraryResultsAutoCompleter->generateProperties(arbitraryResults, structure);

  // Should have computed thermochemical properties
  ASSERT_TRUE(arbitraryResults.has<Property::Thermochemistry>());
  // Insufficient information for calculating Atomic Charges so they should not be there
  ASSERT_FALSE(arbitraryResults.has<Property::AtomicCharges>());

  // Check whether computed Thermochemical properties match orca_test_calc.out
  ASSERT_THAT(arbitraryResults.get<Property::Thermochemistry>().vibrationalComponent.zeroPointVibrationalEnergy,
              DoubleNear(0.02101209, 1e-5));
  ASSERT_THAT(arbitraryResults.get<Property::Thermochemistry>().overall.enthalpy, DoubleNear(-75.79347127, 1e-5));
  ASSERT_THAT(arbitraryResults.get<Property::Thermochemistry>().overall.entropy * 298.15, DoubleNear(0.02210368, 1e-5));
  ASSERT_THAT(arbitraryResults.get<Property::Thermochemistry>().overall.gibbsFreeEnergy, DoubleNear(-75.81557495, 1e-5));
}

TEST_F(AResultsAutoCompleterTest, CanComputeDensityMatrix) {
  // Set required properties
  arbitraryResults.set<Property::SuccessfulCalculation>(true);
  arbitraryResults.set<Property::CoefficientMatrix>(arbitraryCoefficientMatrix);
  arbitraryResults.set<Property::ElectronicOccupation>(arbitraryOccupation);

  arbitraryResultsAutoCompleter = std::make_unique<ResultsAutoCompleter>(randomAtomCollection);

  // Should be able to compute DensityMatrix
  ASSERT_TRUE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::DensityMatrix));
  // Launch calculation
  arbitraryResultsAutoCompleter->generateProperties(arbitraryResults, randomAtomCollection);
  // Should have computed density matrix
  ASSERT_TRUE(arbitraryResults.has<Property::DensityMatrix>());
}

TEST_F(AResultsAutoCompleterTest, CanComputeAtomicChargesAndBondOrderMatrix) {
  // Set required properties
  arbitraryResults.set<Property::SuccessfulCalculation>(true);
  arbitraryResults.set<Property::CoefficientMatrix>(arbitraryCoefficientMatrix); // used for DensityMatrix
  arbitraryResults.set<Property::ElectronicOccupation>(arbitraryOccupation);     // used for DensityMatrix
  arbitraryResults.set<Property::OverlapMatrix>(arbitraryOverlapMatrix);
  arbitraryResults.set<Property::AOtoAtomMapping>(arbitraryAOToAtomMapping);

  arbitraryResultsAutoCompleter = std::make_unique<ResultsAutoCompleter>(randomAtomCollection);

  // Should be able to compute DensityMatrix
  ASSERT_TRUE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::DensityMatrix));
  // Should not be able to compute AtomicCharges and BondOrders right now (depend on DensityMatrix)
  ASSERT_FALSE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::AtomicCharges));
  ASSERT_FALSE(arbitraryResultsAutoCompleter->propertyGeneratable(arbitraryResults, Property::BondOrderMatrix));
  // Launch calculation
  arbitraryResultsAutoCompleter->generateProperties(arbitraryResults, randomAtomCollection);

  // Should have computed density matrix, atomic charges and bond order matrix
  ASSERT_TRUE(arbitraryResults.has<Property::DensityMatrix>());
  ASSERT_TRUE(arbitraryResults.has<Property::AtomicCharges>());
  ASSERT_TRUE(arbitraryResults.has<Property::BondOrderMatrix>());
}

} // namespace Utils
} // namespace Scine
