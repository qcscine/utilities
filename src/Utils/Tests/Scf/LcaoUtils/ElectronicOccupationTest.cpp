/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "gmock/gmock.h"
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace LcaoUtils {
namespace Tests {

class AElectronicOccupation : public Test {
 public:
  ElectronicOccupation occupation;
  int randomN = 134;
  int randomNA = 18;
  int randomNB = 45;
};

TEST_F(AElectronicOccupation, ReturnsCorrectNumberOfElectronsWhenSettingFromNumberOfElectronsForRHF) {
  occupation.fillLowestRestrictedOrbitalsWithElectrons(randomN);

  ASSERT_TRUE(occupation.isFilledUpFromTheBottom());
  ASSERT_THAT(occupation.numberRestrictedElectrons(), Eq(randomN));
  ASSERT_THAT(occupation.numberAlphaElectrons(), Eq(0));
  ASSERT_THAT(occupation.numberBetaElectrons(), Eq(0));
}

TEST_F(AElectronicOccupation, ReturnsCorrectNumberOfElectronsWhenSettingFromNumberOfElectronsForUHF) {
  occupation.fillLowestUnrestrictedOrbitals(randomNA, randomNB);

  ASSERT_TRUE(occupation.isFilledUpFromTheBottom());
  ASSERT_THAT(occupation.numberRestrictedElectrons(), Eq(0));
  ASSERT_THAT(occupation.numberAlphaElectrons(), Eq(randomNA));
  ASSERT_THAT(occupation.numberBetaElectrons(), Eq(randomNB));
}

TEST_F(AElectronicOccupation, ReturnsCorrectNumberOfElectronsWhenOrbitalsSpecified) {
  std::vector<int> orbitals = {0, 2, 3};
  occupation.fillSpecifiedRestrictedOrbitals(orbitals);
  ASSERT_THAT(occupation.numberRestrictedElectrons(), Eq(6));

  std::vector<int> alphaOrbitals = {0, 2, 3};
  std::vector<int> betaOrbitals = {3};
  occupation.fillSpecifiedUnrestrictedOrbitals(alphaOrbitals, betaOrbitals);
  ASSERT_THAT(occupation.numberAlphaElectrons(), Eq(3));
  ASSERT_THAT(occupation.numberBetaElectrons(), Eq(1));
}

TEST_F(AElectronicOccupation, IsOverwrittenWhenReassigned) {
  occupation.fillLowestRestrictedOrbitalsWithElectrons(randomN);
  occupation.fillLowestUnrestrictedOrbitals(randomNA, randomNB);

  ASSERT_THAT(occupation.numberRestrictedElectrons(), Eq(0));

  occupation.fillLowestRestrictedOrbitalsWithElectrons(randomN);

  ASSERT_THAT(occupation.numberAlphaElectrons(), Eq(0));
  ASSERT_THAT(occupation.numberBetaElectrons(), Eq(0));
}

TEST_F(AElectronicOccupation, CanReturnArrayOfOccupiedOrbitalsWhenFilledUpFromBelowForRHF) {
  occupation.fillLowestRestrictedOrbitalsWithElectrons(10);

  auto occupiedOrbitals = occupation.getFilledRestrictedOrbitals();

  decltype(occupiedOrbitals) expected = {0, 1, 2, 3, 4};
  ASSERT_THAT(occupiedOrbitals, Eq(expected));
}

TEST_F(AElectronicOccupation, CanReturnArrayOfOccupiedOrbitalsWhenFilledUpFromBelowForUHF) {
  occupation.fillLowestUnrestrictedOrbitals(3, 1);

  auto occupiedOrbitalsA = occupation.getFilledAlphaOrbitals();
  auto occupiedOrbitalsB = occupation.getFilledBetaOrbitals();

  decltype(occupiedOrbitalsA) expectedA = {0, 1, 2};
  decltype(occupiedOrbitalsB) expectedB = {0};
  ASSERT_THAT(occupiedOrbitalsA, Eq(expectedA));
  ASSERT_THAT(occupiedOrbitalsB, Eq(expectedB));
}

TEST_F(AElectronicOccupation, DetectsWhetherAllConsecutiveWhenOrbitalsAreExplicitlySpecified) {
  occupation.fillSpecifiedRestrictedOrbitals({0, 1, 2, 3, 4});
  ASSERT_TRUE(occupation.isFilledUpFromTheBottom());

  occupation.fillSpecifiedRestrictedOrbitals({0, 1, 2, 4});
  ASSERT_FALSE(occupation.isFilledUpFromTheBottom());

  occupation.fillSpecifiedUnrestrictedOrbitals({0, 1, 2, 3, 4}, {0, 1});
  ASSERT_TRUE(occupation.isFilledUpFromTheBottom());

  occupation.fillSpecifiedUnrestrictedOrbitals({0, 1, 2, 3}, {0, 2});
  ASSERT_FALSE(occupation.isFilledUpFromTheBottom());
}

TEST_F(AElectronicOccupation, CanBeMadeUnrestricted) {
  std::vector<int> occupiedOrbitals = {0, 2, 3, 4};
  occupation.fillSpecifiedRestrictedOrbitals(occupiedOrbitals);

  occupation.makeUnrestricted();

  ASSERT_THAT(occupation.getFilledAlphaOrbitals(), Eq(occupiedOrbitals));
  ASSERT_THAT(occupation.getFilledBetaOrbitals(), Eq(occupiedOrbitals));
}
} // namespace Tests
} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine
