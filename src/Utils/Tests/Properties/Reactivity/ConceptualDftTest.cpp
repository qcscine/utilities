/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Properties/Reactivity/ConceptualDft.h>
#include <gmock/gmock.h>
#include <Eigen/Core>
#include <memory>

using namespace testing;

namespace Scine {
namespace Utils {

/*
 * Reference calculation performed Multiwfn Version 3.7(dev)
 * http://sobereva.com/multiwfn
 * Tian Lu, Feiwu Chen, J. Comput. Chem. 33, 580-592 (2012).
 * ---
 *
 * Note: the E(HOMO) of TCE used for evaluating nucleophilicity index is the value evaluated at B3LYP/6-31G* level
 *
 * Hirshfeld charges, condensed Fukui functions and condensed dual descriptors
 * Units used below are "e" (elementary charge)
 *     Atom     q(N)    q(N+1)   q(N-1)     f-       f+       f0      CDD
 *     1(C )   0.1326  -0.0595   0.2626   0.1301   0.1921   0.1611   0.0621
 *     2(C )   0.1326  -0.0595   0.2624   0.1298   0.1921   0.1610   0.0623
 *     3(H )   0.0419  -0.0547   0.1615   0.1196   0.0966   0.1081  -0.0230
 *     4(O )  -0.1744  -0.3857   0.0751   0.2495   0.2113   0.2304  -0.0382
 *     5(O )  -0.1744  -0.3858   0.0769   0.2513   0.2113   0.2313  -0.0400
 *     6(H )   0.0419  -0.0547   0.1615   0.1197   0.0966   0.1081  -0.0231

 * Condensed local electrophilicity/nucleophilicity index (e*eV)
 *     Atom              Electrophilicity          Nucleophilicity
 *     1(C )                  0.27323                  0.39891
 *     2(C )                  0.27322                  0.39829
 *     3(H )                  0.13733                  0.36675
 *     4(O )                  0.30056                  0.76544
 *     5(O )                  0.30057                  0.77094
 *     6(H )                  0.13733                  0.36707

 * Condensed local softnesses (Hartree*e) and relative electrophilicity/nucleophilicity (dimensionless)
 *     Atom         s-          s+          s0        s+/s-       s-/s+
 *     1(C )      0.3785      0.5591      0.4688      1.4772      0.6770
 *     2(C )      0.3779      0.5591      0.4685      1.4795      0.6759
 *     3(H )      0.3480      0.2810      0.3145      0.8076      1.2383
 *     4(O )      0.7263      0.6151      0.6707      0.8469      1.1808
 *     5(O )      0.7315      0.6151      0.6733      0.8408      1.1893
 *     6(H )      0.3483      0.2810      0.3147      0.8069      1.2393

 * E(N):     -227.402193 Hartree
 * E(N+1):   -227.419914 Hartree
 * E(N-1):   -227.040885 Hartree
 * E_HOMO(N):     -0.222474 Hartree,   -6.0538 eV
 * E_HOMO(N+1):    0.116718 Hartree,    3.1761 eV
 * E_HOMO(N-1):   -0.538957 Hartree,  -14.6658 eV
 * Vertical IP:    0.361308 Hartree,    9.8317 eV
 * Vertical EA:    0.017722 Hartree,    0.4822 eV
 * Mulliken electronegativity:     0.189515 Hartree,    5.1570 eV
 * Chemical potential:            -0.189515 Hartree,   -5.1570 eV
 * Hardness (=fundamental gap):    0.343586 Hartree,    9.3495 eV
 * Softness:    2.910477 Hartree^-1,    0.1070 eV^-1
 * Electrophilicity index:    0.052266 Hartree,    1.4222 eV
 * Nucleophilicity index:     0.112724 Hartree,    3.0674 eV
 */

class AConceptualDftTest : public Test {
 protected:
  // All values in atomic units
  // Input data
  double energy = -227.402193;
  double energyPlus = -227.419914;
  double energyMinus = -227.040885;
  Eigen::VectorXd atomicCharges;
  Eigen::VectorXd atomicChargesPlus;
  Eigen::VectorXd atomicChargesMinus;
  // Multiwfn results
  double chemicalPotential = -0.189515;
  double electronegativity = 0.189515;
  double hardness = 0.343586;
  double softness = 2.910477;
  double electrophilicity = 0.052266;
  Eigen::VectorXd fukuiPlus;
  Eigen::VectorXd fukuiMinus;
  Eigen::VectorXd fukuiRadical;
  Eigen::VectorXd dualDescriptor;

  void SetUp() final {
    // Input data
    atomicCharges.resize(6);
    atomicCharges << 0.1326, 0.1326, 0.0419, -0.1744, -0.1744, 0.0419;
    atomicChargesPlus.resize(6);
    atomicChargesPlus << -0.0595, -0.0595, -0.0547, -0.3857, -0.3858, -0.0547;
    atomicChargesMinus.resize(6);
    atomicChargesMinus << 0.2626, 0.2624, 0.1615, 0.0751, 0.0769, 0.1615;
    // Multiwfn results
    fukuiPlus.resize(6);
    fukuiPlus << 0.1921, 0.1921, 0.0966, 0.2113, 0.2113, 0.0966;
    fukuiMinus.resize(6);
    fukuiMinus << 0.1301, 0.1298, 0.1196, 0.2495, 0.2513, 0.1197;
    fukuiRadical.resize(6);
    fukuiRadical << 0.1611, 0.1610, 0.1081, 0.2304, 0.2313, 0.1081;
    dualDescriptor.resize(6);
    dualDescriptor << 0.0621, 0.0623, -0.0230, -0.0382, -0.0400, -0.0231;
  }
};

TEST_F(AConceptualDftTest, CorrectlyCalculatesChemicalPotential) {
  double calculatedChemicalPotential = ConceptualDft::calculateChemicalPotential(energy, energyPlus, energyMinus);
  ASSERT_THAT(calculatedChemicalPotential, DoubleNear(chemicalPotential, 1e-5));
}
TEST_F(AConceptualDftTest, CorrectlyCalculatesElectronegativity) {
  double calculatedElectronegativity = ConceptualDft::calculateElectronegativity(energy, energyPlus, energyMinus);
  ASSERT_THAT(calculatedElectronegativity, DoubleNear(electronegativity, 1e-5));
}
TEST_F(AConceptualDftTest, CorrectlyCalculatesHardness) {
  double calculatedHardness = ConceptualDft::calculateHardness(energy, energyPlus, energyMinus);
  ASSERT_THAT(calculatedHardness, DoubleNear(hardness, 1e-5));
}
TEST_F(AConceptualDftTest, CorrectlyCalculatesSoftness) {
  double calculatedSoftness = ConceptualDft::calculateSoftness(energy, energyPlus, energyMinus);
  ASSERT_THAT(calculatedSoftness, DoubleNear(softness, 1e-5));
}
TEST_F(AConceptualDftTest, CorrectlyCalculatesElectrophilicity) {
  double calculatedElectrophilicity = ConceptualDft::calculateElectrophilicity(energy, energyPlus, energyMinus);
  ASSERT_THAT(calculatedElectrophilicity, DoubleNear(electrophilicity, 1e-5));
}

TEST_F(AConceptualDftTest, CorrectlyCalculatesFukuiPlus) {
  Eigen::VectorXd calculatedFukuiPlus =
      ConceptualDft::calculateFukuiPlus(atomicCharges, atomicChargesPlus, atomicChargesMinus);
  ASSERT_FALSE(((calculatedFukuiPlus - fukuiPlus).array().abs() > 1e-4).any());
}

TEST_F(AConceptualDftTest, CorrectlyCalculatesFukuiMinus) {
  Eigen::VectorXd calculatedFukuiMinus =
      ConceptualDft::calculateFukuiMinus(atomicCharges, atomicChargesPlus, atomicChargesMinus);
  ASSERT_FALSE(((calculatedFukuiMinus - fukuiMinus).array().abs() > 1e-4).any());
}

TEST_F(AConceptualDftTest, CorrectlyCalculatesFukuiRadical) {
  Eigen::VectorXd calculatedFukuiRadical =
      ConceptualDft::calculateFukuiRadical(atomicCharges, atomicChargesPlus, atomicChargesMinus);
  ASSERT_FALSE(((calculatedFukuiRadical - fukuiRadical).array().abs() > 1e-4).any());
}

TEST_F(AConceptualDftTest, CorrectlyCalculatesDualDescriptor) {
  Eigen::VectorXd calculatedDualDescriptor =
      ConceptualDft::calculateDualDescriptor(atomicCharges, atomicChargesPlus, atomicChargesMinus);
  ASSERT_FALSE(((calculatedDualDescriptor - dualDescriptor).array().abs() > 2e-4).any());
}

TEST_F(AConceptualDftTest, CorrectlyFillsGlobalContainer) {
  auto globalContainer = ConceptualDft::calculateGlobal(energy, energyPlus, energyMinus);
  ASSERT_THAT(globalContainer.chemicalPotential, DoubleNear(chemicalPotential, 1e-5));
  ASSERT_THAT(globalContainer.electronegativity, DoubleNear(electronegativity, 1e-5));
  ASSERT_THAT(globalContainer.hardness, DoubleNear(hardness, 1e-5));
  ASSERT_THAT(globalContainer.softness, DoubleNear(softness, 1e-5));
  ASSERT_THAT(globalContainer.electrophilicity, DoubleNear(electrophilicity, 1e-5));
}

TEST_F(AConceptualDftTest, CorrectlyFillsLocalContainer) {
  auto localContainer = ConceptualDft::calculateLocal(atomicCharges, atomicChargesPlus, atomicChargesMinus);
  ASSERT_FALSE(((localContainer.fukuiPlus - fukuiPlus).array().abs() > 1e-4).any());
  ASSERT_FALSE(((localContainer.fukuiMinus - fukuiMinus).array().abs() > 1e-4).any());
  ASSERT_FALSE(((localContainer.fukuiRadical - fukuiRadical).array().abs() > 1e-4).any());
  ASSERT_FALSE(((localContainer.dualDescriptor - dualDescriptor).array().abs() > 2e-4).any());
}

TEST_F(AConceptualDftTest, CorrectlyFillsOverallContainer) {
  auto overallContainer =
      ConceptualDft::calculate(energy, atomicCharges, energyPlus, atomicChargesPlus, energyMinus, atomicChargesMinus);
  ASSERT_THAT(overallContainer.global.chemicalPotential, DoubleNear(chemicalPotential, 1e-5));
  ASSERT_THAT(overallContainer.global.electronegativity, DoubleNear(electronegativity, 1e-5));
  ASSERT_THAT(overallContainer.global.hardness, DoubleNear(hardness, 1e-5));
  ASSERT_THAT(overallContainer.global.softness, DoubleNear(softness, 1e-5));
  ASSERT_THAT(overallContainer.global.electrophilicity, DoubleNear(electrophilicity, 1e-5));
  ASSERT_FALSE(((overallContainer.local.fukuiPlus - fukuiPlus).array().abs() > 1e-4).any());
  ASSERT_FALSE(((overallContainer.local.fukuiMinus - fukuiMinus).array().abs() > 1e-4).any());
  ASSERT_FALSE(((overallContainer.local.fukuiRadical - fukuiRadical).array().abs() > 1e-4).any());
  ASSERT_FALSE(((overallContainer.local.dualDescriptor - dualDescriptor).array().abs() > 2e-4).any());
}

} // namespace Utils
} // namespace Scine
