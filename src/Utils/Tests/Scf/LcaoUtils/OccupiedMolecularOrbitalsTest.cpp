/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/OccupiedMolecularOrbitals.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {
class AOccupiedMolecularOrbitals : public Test {
 public:
  MolecularOrbitals randomRhfOrbitals;
  MolecularOrbitals randomUhfOrbitals;
  const int nAOs = 8;

  void SetUp() override {
    Eigen::MatrixXd randomMatrix1 = Eigen::MatrixXd::Random(nAOs, nAOs);
    Eigen::MatrixXd randomMatrix2 = Eigen::MatrixXd::Random(nAOs, nAOs);
    Eigen::MatrixXd randomMatrix3 = Eigen::MatrixXd::Random(nAOs, nAOs);

    randomRhfOrbitals = MolecularOrbitals::createFromRestrictedCoefficients(std::move(randomMatrix1));
    randomUhfOrbitals = MolecularOrbitals::createFromUnrestrictedCoefficients(randomMatrix2, randomMatrix3);
  }
};

TEST_F(AOccupiedMolecularOrbitals, IsEquivalentToMolecularOrbitalsIfAllOrbitalsAreFilled) {
  LcaoUtils::ElectronicOccupation fullRhfOccupation;
  LcaoUtils::ElectronicOccupation fullUhfOccupation;
  fullRhfOccupation.fillLowestRestrictedOrbitalsWithElectrons(2 * nAOs);
  fullUhfOccupation.fillLowestUnrestrictedOrbitals(nAOs, nAOs);

  OccupiedMolecularOrbitals occRhf{randomRhfOrbitals, fullRhfOccupation};
  OccupiedMolecularOrbitals occUhf{randomUhfOrbitals, fullUhfOccupation};

  ASSERT_TRUE(occRhf.restrictedMatrix().isApprox(randomRhfOrbitals.restrictedMatrix()));
  ASSERT_TRUE(occUhf.alphaMatrix().isApprox(randomUhfOrbitals.alphaMatrix()));
  ASSERT_TRUE(occUhf.betaMatrix().isApprox(randomUhfOrbitals.betaMatrix()));
}

TEST_F(AOccupiedMolecularOrbitals, ConstructedCorrectlyForRandomRhfOccupation) {
  LcaoUtils::ElectronicOccupation rhfOccupation;
  rhfOccupation.fillSpecifiedRestrictedOrbitals({0, 3, 5});

  OccupiedMolecularOrbitals occRhf{randomRhfOrbitals, rhfOccupation};

  ASSERT_TRUE(occRhf.isRestricted());
  ASSERT_TRUE(occRhf.restrictedMatrix().col(0).isApprox(randomRhfOrbitals.restrictedMatrix().col(0)));
  ASSERT_TRUE(occRhf.restrictedMatrix().col(1).isApprox(randomRhfOrbitals.restrictedMatrix().col(3)));
  ASSERT_TRUE(occRhf.restrictedMatrix().col(2).isApprox(randomRhfOrbitals.restrictedMatrix().col(5)));
}

TEST_F(AOccupiedMolecularOrbitals, ConstructedCorrectlyForRandomUhfOccupation) {
  LcaoUtils::ElectronicOccupation uhfOccupation;
  uhfOccupation.fillSpecifiedUnrestrictedOrbitals({0, 3, 5, 6}, {1, 3});

  OccupiedMolecularOrbitals occUhf{randomUhfOrbitals, uhfOccupation};

  ASSERT_TRUE(occUhf.isUnrestricted());
  ASSERT_TRUE(occUhf.alphaMatrix().col(0).isApprox(randomUhfOrbitals.alphaMatrix().col(0)));
  ASSERT_TRUE(occUhf.alphaMatrix().col(1).isApprox(randomUhfOrbitals.alphaMatrix().col(3)));
  ASSERT_TRUE(occUhf.alphaMatrix().col(2).isApprox(randomUhfOrbitals.alphaMatrix().col(5)));
  ASSERT_TRUE(occUhf.alphaMatrix().col(3).isApprox(randomUhfOrbitals.alphaMatrix().col(6)));
  ASSERT_TRUE(occUhf.betaMatrix().col(0).isApprox(randomUhfOrbitals.betaMatrix().col(1)));
  ASSERT_TRUE(occUhf.betaMatrix().col(1).isApprox(randomUhfOrbitals.betaMatrix().col(3)));
}
} // namespace Tests
} // namespace Utils
} // namespace Scine
