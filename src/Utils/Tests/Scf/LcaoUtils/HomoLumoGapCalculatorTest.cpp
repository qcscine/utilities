/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Scf/LcaoUtils/HomoLumoGapCalculator.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {
class AHomoLumoGapCalculator : public Test {
 public:
  SingleParticleEnergies energies;
  LcaoUtils::ElectronicOccupation occupation;

  void SetUp() override {
  }
};

TEST_F(AHomoLumoGapCalculator, CalculatesHomoLumoGapForRestrictedWithEvenElectrons) {
  Eigen::VectorXd values(6);
  values << -2.2, 0, 3.3, 9.9, 30, 40.2;

  energies.setRestricted(values);
  occupation.fillLowestRestrictedOrbitalsWithElectrons(6);

  double homoLumoGap = LcaoUtils::HomoLumoGapCalculator::calculate(energies, occupation);

  ASSERT_THAT(homoLumoGap, DoubleEq(9.9 - 3.3));
}

TEST_F(AHomoLumoGapCalculator, CalculatesHomoLumoGapForRestrictedWithOddElectrons) {
  Eigen::VectorXd values(6);
  values << -2.2, 0, 3.3, 9.9, 30, 40.2;

  energies.setRestricted(values);
  occupation.fillLowestRestrictedOrbitalsWithElectrons(7);

  double homoLumoGap = LcaoUtils::HomoLumoGapCalculator::calculate(energies, occupation);

  ASSERT_THAT(homoLumoGap, DoubleEq(30 - 9.9));
}

TEST_F(AHomoLumoGapCalculator, CalculatesHomoLumoGapForUnrestricted) {
  Eigen::VectorXd valuesAlpha(5), valuesBeta(5);
  valuesAlpha << -2.2, 0, 3.3, 9.9, 10;
  valuesBeta << -3.2, -1, 10.3, 12.9, 20;

  energies.setUnrestricted(valuesAlpha, valuesBeta);
  occupation.fillLowestUnrestrictedOrbitals(3, 2);

  double homoLumoGap = LcaoUtils::HomoLumoGapCalculator::calculate(energies, occupation);

  ASSERT_THAT(homoLumoGap, DoubleEq(9.9 - 3.3));
}

TEST_F(AHomoLumoGapCalculator, CalculatesHomoLumoGapForUnrestrictedWithBetaPolarization) {
  Eigen::VectorXd valuesAlpha(5), valuesBeta(5);
  valuesAlpha << -3.2, -1, 10.3, 12.9, 20;
  valuesBeta << -2.2, 0, 3.3, 9.9, 10;

  energies.setUnrestricted(valuesAlpha, valuesBeta);
  occupation.fillLowestUnrestrictedOrbitals(2, 3);

  double homoLumoGap = LcaoUtils::HomoLumoGapCalculator::calculate(energies, occupation);

  ASSERT_THAT(homoLumoGap, DoubleEq(9.9 - 3.3));
}

TEST_F(AHomoLumoGapCalculator, CalculatesHomoLumoGapForUnrestrictedAndHOMOLUMOWithDifferentPolarity) {
  Eigen::VectorXd valuesAlpha(5), valuesBeta(5);
  valuesAlpha << -2.2, 0, 3.3, 10.3, 10;
  valuesBeta << -3.2, -1, 9.9, 12.9, 20;

  energies.setUnrestricted(valuesAlpha, valuesBeta);
  occupation.fillLowestUnrestrictedOrbitals(3, 2);

  double homoLumoGap = LcaoUtils::HomoLumoGapCalculator::calculate(energies, occupation);

  ASSERT_THAT(homoLumoGap, DoubleEq(9.9 - 3.3));
}

TEST_F(AHomoLumoGapCalculator, ThrowsWhenThereAreNoElectrons) {
  Eigen::VectorXd values(6);
  values << -2.2, 0, 3.3, 9.9, 30, 40.2;
  energies.setRestricted(values);
  occupation.fillLowestRestrictedOrbitalsWithElectrons(0);
  ASSERT_THROW(LcaoUtils::HomoLumoGapCalculator::calculate(energies, occupation), LcaoUtils::HomoLumoGapException);

  Eigen::VectorXd valuesAlpha(5), valuesBeta(5);
  valuesAlpha << -2.2, 0, 3.3, 9.9, 10;
  valuesBeta << -3.2, -1, 10.3, 12.9, 20;
  energies.setUnrestricted(valuesAlpha, valuesBeta);
  occupation.fillLowestUnrestrictedOrbitals(0, 0);
  ASSERT_THROW(LcaoUtils::HomoLumoGapCalculator::calculate(energies, occupation), LcaoUtils::HomoLumoGapException);
}

TEST_F(AHomoLumoGapCalculator, ThrowsWhenThereAreNoEmptyOrbitals) {
  Eigen::VectorXd values(6);
  values << -2.2, 0, 3.3, 9.9, 30, 40.2;
  energies.setRestricted(values);
  occupation.fillLowestRestrictedOrbitalsWithElectrons(12);
  ASSERT_THROW(LcaoUtils::HomoLumoGapCalculator::calculate(energies, occupation), LcaoUtils::HomoLumoGapException);

  Eigen::VectorXd valuesAlpha(5), valuesBeta(5);
  valuesAlpha << -2.2, 0, 3.3, 9.9, 10;
  valuesBeta << -3.2, -1, 10.3, 12.9, 20;
  energies.setUnrestricted(valuesAlpha, valuesBeta);
  occupation.fillLowestUnrestrictedOrbitals(5, 5);
  ASSERT_THROW(LcaoUtils::HomoLumoGapCalculator::calculate(energies, occupation), LcaoUtils::HomoLumoGapException);
}
} // namespace Tests
} // namespace Utils
} // namespace Scine
