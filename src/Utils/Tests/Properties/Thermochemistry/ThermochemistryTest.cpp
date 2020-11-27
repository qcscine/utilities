/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/GeometricDerivatives/NormalModeAnalysis.h>
#include <Utils/Geometry.h>
#include <Utils/Properties/Thermochemistry/ThermochemistryCalculator.h>
#include <Utils/Typenames.h>
#include <gmock/gmock.h>
#include <memory>

using namespace testing;

namespace Scine {
namespace Utils {

/*
 * Reference calculation performed with MOPAC2016.
 * PM6 PRECISE
 *
 *
 * C     0.000000  1   0.000000  1   0.000000 1
 * O     0.000000  1   0.000000  1   1.212200 1
 * H     0.937197  1   0.000000  1  -0.584262 1
 * H    -0.937197  1   0.000000  1  -0.584262 1
 *
 * PM6 RHF THERMO OLDGEO PRECISE
 *
 * and
 *
 * PM6 PRECISE
 *
 *
 * H     0  1   0.000000  0  0 0
 * F     1  1   0.000000  0  0 0
 *
 * PM6 RHF THERMO OLDGEO PRECISE
 */

class AThermochemistryTest : public Test {
 public:
  NormalModesContainer arbitraryNormalModes;
  Geometry::PrincipalMomentsOfInertia arbitraryPMI;
  std::unique_ptr<ThermochemistryCalculator> arbitraryTCCalculator;
  ElementTypeCollection arbitraryElements;
  int arbitraryMultiplicity = 1;
  double arbitraryEnergy = 1.0;

 protected:
  void SetUp() final {
    arbitraryElements = {ElementType::C, ElementType::O, ElementType::H, ElementType::H};
    NormalMode m1(1101.75, DisplacementCollection::Random(4, 3));
    arbitraryNormalModes.add(m1);
    NormalMode m2(1157.94, DisplacementCollection::Random(4, 3));
    arbitraryNormalModes.add(m2);
    NormalMode m3(1349.03, DisplacementCollection::Random(4, 3));
    arbitraryNormalModes.add(m3);
    NormalMode m4(1791.24, DisplacementCollection::Random(4, 3));
    arbitraryNormalModes.add(m4);
    NormalMode m5(2614.79, DisplacementCollection::Random(4, 3));
    arbitraryNormalModes.add(m5);
    NormalMode m6(2664.54, DisplacementCollection::Random(4, 3));
    arbitraryNormalModes.add(m6);

    Eigen::Vector3d eigenValues(2.8969, 21.7672, 24.6640);
    // Convert from 1e-40 g cm^2 to amu*bohr^2
    eigenValues *= 1e-47 * Constants::u_per_kg * std::pow(Constants::bohr_per_meter, 2);
    arbitraryPMI.eigenvalues = eigenValues;
    arbitraryPMI.eigenvectors = Eigen::Matrix3d::Random();
  }
};

TEST_F(AThermochemistryTest, CanConstructCalculator) {
  arbitraryTCCalculator = std::make_unique<ThermochemistryCalculator>(arbitraryNormalModes, arbitraryPMI, arbitraryElements,
                                                                      arbitraryMultiplicity, arbitraryEnergy);
}

TEST_F(AThermochemistryTest, CorrectlyCalculatesZPVE) {
  arbitraryTCCalculator = std::make_unique<ThermochemistryCalculator>(arbitraryNormalModes, arbitraryPMI, arbitraryElements,
                                                                      arbitraryMultiplicity, arbitraryEnergy);
  arbitraryTCCalculator->setZPVEInclusion(ZPVEInclusion::alreadyIncluded);
  arbitraryTCCalculator->setTemperature(298.15);
  auto container = arbitraryTCCalculator->calculate();
  ASSERT_THAT(container.overall.symmetryNumber, Eq(1));
  ASSERT_THAT(container.vibrationalComponent.zeroPointVibrationalEnergy * Constants::kCalPerMol_per_hartree,
              DoubleNear(15.267, 1e-3));
}

TEST_F(AThermochemistryTest, CorrectlyCalculatesZPVEForLinearMolecule) {
  NormalModesContainer HFNormalModes;
  NormalMode m1(3968.7, DisplacementCollection::Random(1, 3));
  HFNormalModes.add(m1);
  arbitraryTCCalculator = std::make_unique<ThermochemistryCalculator>(HFNormalModes, arbitraryPMI, arbitraryElements,
                                                                      arbitraryMultiplicity, arbitraryEnergy);
  arbitraryTCCalculator->setZPVEInclusion(ZPVEInclusion::alreadyIncluded);
  arbitraryTCCalculator->setTemperature(298.15);

  auto container = arbitraryTCCalculator->calculate();
  ASSERT_THAT(container.vibrationalComponent.zeroPointVibrationalEnergy * Constants::kCalPerMol_per_hartree,
              DoubleNear(5.674, 1e-3));
}

TEST_F(AThermochemistryTest, CanCalculateVibrationalComponent) {
  arbitraryTCCalculator = std::make_unique<ThermochemistryCalculator>(arbitraryNormalModes, arbitraryPMI, arbitraryElements,
                                                                      arbitraryMultiplicity, arbitraryEnergy);
  arbitraryTCCalculator->setZPVEInclusion(ZPVEInclusion::alreadyIncluded);
  arbitraryTCCalculator->setTemperature(298.00);
  auto container = arbitraryTCCalculator->calculate();

  EXPECT_THAT(container.vibrationalComponent.heatCapacityP * Constants::kCalPerMol_per_hartree,
              DoubleNear(0.6650 / 1000, 1e-6));
  EXPECT_THAT(container.vibrationalComponent.entropy * Constants::kCalPerMol_per_hartree, DoubleNear(0.1365 / 1000, 1e-6));
  EXPECT_THAT(container.vibrationalComponent.enthalpy * Constants::kCalPerMol_per_hartree, DoubleNear(0.0345759, 1e-6));
}

TEST_F(AThermochemistryTest, CanCalculateRotationalComponentHF) {
  NormalModesContainer HFNormalModes;
  NormalMode m1(3968.7, DisplacementCollection::Random(1, 3));
  HFNormalModes.add(m1);

  ElementTypeCollection elements = {ElementType::H, ElementType::F};

  Eigen::Vector3d eigenValues(0.00000000, 1.4818, 1.4818);
  // Convert from 1e-40 g cm^2 to amu*bohr^2
  eigenValues *= 1e-47 * Constants::u_per_kg * std::pow(Constants::bohr_per_meter, 2);
  Eigen::Matrix3d eigenVectors = Eigen::Matrix3d::Random();
  Geometry::PrincipalMomentsOfInertia pmi;
  pmi.eigenvalues = eigenValues;
  pmi.eigenvectors = eigenVectors;

  arbitraryTCCalculator = std::make_unique<ThermochemistryCalculator>(HFNormalModes, pmi, elements, 1, arbitraryEnergy);
  arbitraryTCCalculator->setZPVEInclusion(ZPVEInclusion::alreadyIncluded);
  arbitraryTCCalculator->setTemperature(298.00);
  auto container = arbitraryTCCalculator->calculate();
  ASSERT_THAT(container.rotationalComponent.enthalpy * Constants::kCalPerMol_per_hartree, DoubleNear(592.1875 / 1000, 1e-5));
  ASSERT_THAT(container.rotationalComponent.heatCapacityP * Constants::kCalPerMol_per_hartree, DoubleNear(1.9872 / 1000, 1e-5));
  ASSERT_THAT(container.rotationalComponent.entropy * Constants::kCalPerMol_per_hartree, DoubleNear(6.7458 / 1000, 1e-5));
}

TEST_F(AThermochemistryTest, CanCalculateRotationalComponentOfAsymmetricMolecule) {
  Eigen::Matrix3d eigenVectors = Eigen::Matrix3d::Random();

  arbitraryTCCalculator = std::make_unique<ThermochemistryCalculator>(arbitraryNormalModes, arbitraryPMI, arbitraryElements,
                                                                      arbitraryMultiplicity, arbitraryEnergy);
  arbitraryTCCalculator->setZPVEInclusion(ZPVEInclusion::alreadyIncluded);
  arbitraryTCCalculator->setTemperature(298.00);
  // C2v symmetry
  arbitraryTCCalculator->setMolecularSymmetryNumber(2);
  auto container = arbitraryTCCalculator->calculate();
  ASSERT_THAT(container.overall.symmetryNumber, Eq(2));
  ASSERT_THAT(container.rotationalComponent.enthalpy * Constants::kCalPerMol_per_hartree, DoubleNear(888.2813 / 1000, 1e-5));
  ASSERT_THAT(container.rotationalComponent.heatCapacityP * Constants::kCalPerMol_per_hartree, DoubleNear(2.9808 / 1000, 1e-5));
  ASSERT_THAT(container.rotationalComponent.entropy * Constants::kCalPerMol_per_hartree, DoubleNear(16.0088 / 1000, 1e-5));
}

TEST_F(AThermochemistryTest, CanCalculateTranslationalComponentOfAsymmetricMolecule) {
  arbitraryTCCalculator = std::make_unique<ThermochemistryCalculator>(arbitraryNormalModes, arbitraryPMI, arbitraryElements,
                                                                      arbitraryMultiplicity, arbitraryEnergy);
  arbitraryTCCalculator->setZPVEInclusion(ZPVEInclusion::alreadyIncluded);
  arbitraryTCCalculator->setTemperature(298.00);
  // C2v symmetry
  arbitraryTCCalculator->setMolecularSymmetryNumber(2);
  auto container = arbitraryTCCalculator->calculate();
  ASSERT_THAT(container.translationalComponent.enthalpy * Constants::kCalPerMol_per_hartree,
              DoubleNear(1480.4688 / 1000, 1e-5));
  ASSERT_THAT(container.translationalComponent.heatCapacityP * Constants::kCalPerMol_per_hartree,
              DoubleNear(4.9680 / 1000, 1e-5));
  ASSERT_THAT(container.translationalComponent.entropy * Constants::kCalPerMol_per_hartree, DoubleNear(36.1295 / 1000, 1e-5));
}

TEST_F(AThermochemistryTest, CanCalculateOverallComponent) {
  arbitraryTCCalculator = std::make_unique<ThermochemistryCalculator>(arbitraryNormalModes, arbitraryPMI, arbitraryElements,
                                                                      arbitraryMultiplicity, arbitraryEnergy);
  arbitraryTCCalculator->setZPVEInclusion(ZPVEInclusion::alreadyIncluded);
  arbitraryTCCalculator->setTemperature(298.00);
  // C2v symmetry
  arbitraryTCCalculator->setMolecularSymmetryNumber(2);
  auto container = arbitraryTCCalculator->calculate();
  ASSERT_THAT(container.overall.enthalpy * Constants::kCalPerMol_per_hartree,
              DoubleNear(2403.3261 / 1000 + arbitraryEnergy * Constants::kCalPerMol_per_hartree, 3e-5));
  ASSERT_THAT(container.overall.heatCapacityP * Constants::kCalPerMol_per_hartree, DoubleNear(8.6138 / 1000, 3e-5));
  ASSERT_THAT(container.overall.entropy * Constants::kCalPerMol_per_hartree, DoubleNear(52.2772 / 1000, 3e-5));
}
} // namespace Utils
} // namespace Scine
