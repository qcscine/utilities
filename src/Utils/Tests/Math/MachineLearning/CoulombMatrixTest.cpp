/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Math/MachineLearning/ChemicalRepresentations/CoulombMatrix.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
using namespace Utils;
using namespace MachineLearning;
namespace Tests {

/**
 * @class ACoulombMatrixTest CoulombMatrixTest.cpp
 * @class Tests concerning the construction of the Coulomb matrix for an arbitrary molecular system.
 * @test
 */
class ACoulombMatrixTest : public Test {};

TEST_F(ACoulombMatrixTest, IsCorrectlyConstructedForCarbonMonoxide) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000    0.0000000000\n"
                       "O     2.1167100000    0.0000000000    0.0000000000\n");
  auto atomCollection = Utils::XyzStreamHandler::read(ss);

  CoulombMatrix cm(atomCollection);

  auto matrix = cm.getMatrix();
  auto features = cm.getFeatures();

  ASSERT_THAT(features.size(), Eq(3));
  ASSERT_THAT(matrix.size(), Eq(4));

  ASSERT_THAT(matrix(0, 0), DoubleNear(36.8581052, 1e-5));
  ASSERT_THAT(matrix(1, 1), DoubleNear(73.5166947, 1e-5));
  ASSERT_THAT(std::abs(matrix(0, 1) - matrix(1, 0)), DoubleNear(0.0, 1e-5));
  ASSERT_THAT(matrix(0, 1), DoubleNear(12.0, 1e-5));

  ASSERT_THAT(std::abs(matrix(0, 0) - features(0)), DoubleNear(0.0, 1e-5));
  ASSERT_THAT(std::abs(matrix(0, 1) - features(1)), DoubleNear(0.0, 1e-5));
  ASSERT_THAT(std::abs(matrix(1, 1) - features(2)), DoubleNear(0.0, 1e-5));
}

TEST_F(ACoulombMatrixTest, IsCorrectlyConstructedForAlanine) {
  std::stringstream ss("13\n\n"
                       "H   2.14925093178292      0.31949915808954      1.61355371175786\n"
                       "N   1.77782099345575     -0.09932466221128      0.75348143685936\n"
                       "H   2.02850269645183      0.52789790022383     -0.02373591487857\n"
                       "C   0.32747815689941     -0.17807747546200      0.82240049454924\n"
                       "H   0.02257867916421     -1.11256136944157      1.34643074925508\n"
                       "C   -0.37493862961075      1.00578398018676      1.52688386911241\n"
                       "H   -0.01139894586936      1.09147002335672      2.57087134564999\n"
                       "H   -0.14537987535257      1.95493993856176      1.00009318984220\n"
                       "H   -1.47546552849672      0.87180158210824      1.55612322292084\n"
                       "C   -0.21023562319213     -0.26845307119646     -0.60072077306538\n"
                       "O   0.35218728197586      0.17996478649275     -1.58208848064206\n"
                       "O   -1.43081282327720     -0.85313228887232     -0.64945795029584\n"
                       "H   -1.70918731393126     -0.80500850183596     -1.59033490106512\n");
  auto atomCollection = XyzStreamHandler::read(ss);

  CoulombMatrix cm(atomCollection);

  auto matrix = cm.getMatrix();
  auto features = cm.getFeatures();
  auto nAtoms = atomCollection.size();

  ASSERT_THAT(features.size(), Eq(nAtoms * (nAtoms + 1) / 2));
  ASSERT_THAT(matrix.rows(), Eq(nAtoms));
  ASSERT_THAT(matrix.cols(), Eq(nAtoms));

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      ASSERT_THAT(matrix(i, j), DoubleNear(matrix(j, i), 1e-6));
    }
  }

  ASSERT_THAT(matrix(2, 2), DoubleNear(0.5, 1e-3));
  ASSERT_THAT(matrix(1, 10), DoubleNear(10.774, 1e-3));
}

} // namespace Tests
} // namespace Scine
