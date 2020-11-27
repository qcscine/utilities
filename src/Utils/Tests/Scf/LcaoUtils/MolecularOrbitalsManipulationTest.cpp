/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/Scf/LcaoUtils/MolecularOrbitalsManipulation.h>
#include <gmock/gmock.h>
#include <cmath>

using namespace testing;
namespace Scine {
namespace Utils {
namespace LcaoUtils {
namespace Tests {

class AMolecularOrbitalsManipulation : public Test {
 public:
  const int nOrbitals = 10;
  MolecularOrbitals randomRestrictedMOs;
  MolecularOrbitals randomUnrestrictedMOs;

 protected:
  void SetUp() override {
    Eigen::MatrixXd random1 = Eigen::MatrixXd::Random(nOrbitals, nOrbitals);
    Eigen::MatrixXd random2 = Eigen::MatrixXd::Random(nOrbitals, nOrbitals);
    Eigen::MatrixXd random3 = Eigen::MatrixXd::Random(nOrbitals, nOrbitals);

    randomRestrictedMOs = MolecularOrbitals::createFromRestrictedCoefficients(std::move(random1));
    randomUnrestrictedMOs = MolecularOrbitals::createFromUnrestrictedCoefficients(std::move(random2), std::move(random3));
  }

 private:
};

TEST_F(AMolecularOrbitalsManipulation, SwapsRestrictedOrbitalsCorrectly) {
  MolecularOrbitals expected = randomRestrictedMOs;
  std::vector<MolecularOrbitalsManipulation::Swap> swaps = {{4, 1}, {5, 8}};

  expected.restrictedMatrix().col(4) = randomRestrictedMOs.restrictedMatrix().col(1);
  expected.restrictedMatrix().col(1) = randomRestrictedMOs.restrictedMatrix().col(4);
  expected.restrictedMatrix().col(5) = randomRestrictedMOs.restrictedMatrix().col(8);
  expected.restrictedMatrix().col(8) = randomRestrictedMOs.restrictedMatrix().col(5);

  auto newOrbitals = MolecularOrbitalsManipulation::createRestrictedWithSwaps(randomRestrictedMOs, swaps);

  ASSERT_TRUE(newOrbitals.restrictedMatrix().isApprox(expected.restrictedMatrix()));
}

TEST_F(AMolecularOrbitalsManipulation, SwapsUnrestrictedOrbitalsCorrectly) {
  MolecularOrbitals expected = randomUnrestrictedMOs;
  std::vector<MolecularOrbitalsManipulation::Swap> alphaSwaps = {{4, 1}, {5, 8}};
  std::vector<MolecularOrbitalsManipulation::Swap> betaSwaps = {{0, 2}};

  expected.alphaMatrix().col(4) = randomUnrestrictedMOs.alphaMatrix().col(1);
  expected.alphaMatrix().col(1) = randomUnrestrictedMOs.alphaMatrix().col(4);
  expected.alphaMatrix().col(5) = randomUnrestrictedMOs.alphaMatrix().col(8);
  expected.alphaMatrix().col(8) = randomUnrestrictedMOs.alphaMatrix().col(5);
  expected.betaMatrix().col(0) = randomUnrestrictedMOs.betaMatrix().col(2);
  expected.betaMatrix().col(2) = randomUnrestrictedMOs.betaMatrix().col(0);

  auto newOrbitals = MolecularOrbitalsManipulation::createUnrestrictedWithSwaps(randomUnrestrictedMOs, alphaSwaps, betaSwaps);

  ASSERT_TRUE(newOrbitals.alphaMatrix().isApprox(expected.alphaMatrix()));
  ASSERT_TRUE(newOrbitals.betaMatrix().isApprox(expected.betaMatrix()));
}

TEST_F(AMolecularOrbitalsManipulation, MixesRestrictedOrbitalsCorrectly) {
  MolecularOrbitals expected = randomRestrictedMOs;
  std::vector<MolecularOrbitalsManipulation::Mix> mixes = {{4, 1, 0.1234}, {5, 8, 0.854}};

  expected.restrictedMatrix().col(4) = std::cos(0.1234) * randomRestrictedMOs.restrictedMatrix().col(4) +
                                       std::sin(0.1234) * randomRestrictedMOs.restrictedMatrix().col(1);
  expected.restrictedMatrix().col(1) = std::cos(0.1234) * randomRestrictedMOs.restrictedMatrix().col(1) -
                                       std::sin(0.1234) * randomRestrictedMOs.restrictedMatrix().col(4);
  expected.restrictedMatrix().col(5) = std::cos(0.8540) * randomRestrictedMOs.restrictedMatrix().col(5) +
                                       std::sin(0.8540) * randomRestrictedMOs.restrictedMatrix().col(8);
  expected.restrictedMatrix().col(8) = std::cos(0.8540) * randomRestrictedMOs.restrictedMatrix().col(8) -
                                       std::sin(0.8540) * randomRestrictedMOs.restrictedMatrix().col(5);

  auto newOrbitals = MolecularOrbitalsManipulation::createRestrictedWithMixes(randomRestrictedMOs, mixes);

  ASSERT_TRUE(newOrbitals.restrictedMatrix().isApprox(expected.restrictedMatrix()));
}

TEST_F(AMolecularOrbitalsManipulation, MixesUnrestrictedOrbitalsCorrectly) {
  MolecularOrbitals expected = randomUnrestrictedMOs;
  std::vector<MolecularOrbitalsManipulation::Mix> alphaMixes = {{4, 1, 0.1234}, {5, 8, 0.854}};
  std::vector<MolecularOrbitalsManipulation::Mix> betaMixes = {{0, 2, 0.099}};

  expected.alphaMatrix().col(4) = std::cos(0.1234) * randomUnrestrictedMOs.alphaMatrix().col(4) +
                                  std::sin(0.1234) * randomUnrestrictedMOs.alphaMatrix().col(1);
  expected.alphaMatrix().col(1) = std::cos(0.1234) * randomUnrestrictedMOs.alphaMatrix().col(1) -
                                  std::sin(0.1234) * randomUnrestrictedMOs.alphaMatrix().col(4);
  expected.alphaMatrix().col(5) = std::cos(0.854) * randomUnrestrictedMOs.alphaMatrix().col(5) +
                                  std::sin(0.854) * randomUnrestrictedMOs.alphaMatrix().col(8);
  expected.alphaMatrix().col(8) = std::cos(0.854) * randomUnrestrictedMOs.alphaMatrix().col(8) -
                                  std::sin(0.854) * randomUnrestrictedMOs.alphaMatrix().col(5);
  expected.betaMatrix().col(0) = std::cos(0.099) * randomUnrestrictedMOs.betaMatrix().col(0) +
                                 std::sin(0.099) * randomUnrestrictedMOs.betaMatrix().col(2);
  expected.betaMatrix().col(2) = std::cos(0.099) * randomUnrestrictedMOs.betaMatrix().col(2) -
                                 std::sin(0.099) * randomUnrestrictedMOs.betaMatrix().col(0);

  auto newOrbitals = MolecularOrbitalsManipulation::createUnrestrictedWithMixes(randomUnrestrictedMOs, alphaMixes, betaMixes);

  ASSERT_TRUE(newOrbitals.alphaMatrix().isApprox(expected.alphaMatrix()));
  ASSERT_TRUE(newOrbitals.betaMatrix().isApprox(expected.betaMatrix()));
}
} // namespace Tests
} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine
