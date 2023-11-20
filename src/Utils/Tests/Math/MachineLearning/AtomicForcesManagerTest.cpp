/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/MSVCCompatibility.h"
#include <Utils/Dftd3/Dftd3.h>
#include <Utils/Geometry/GeometryUtilities.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Math/MachineLearning/ChemicalRepresentations/AtomicForcesManager.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
using namespace Utils;
using namespace MachineLearning;
namespace Tests {

/**
 * @class AnAtomicForcesManagerTest AtomicForcesManagerTest.cpp
 * @class Tests the functionalities of the AtomicForcesManager class.
 * @test
 */
class AnAtomicForcesManagerTest : public Test {};

TEST_F(AnAtomicForcesManagerTest, InternalForceRepresentationsAreRotationallyInvariant) {
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
  auto structure = XyzStreamHandler::read(ss);

  AtomicForcesManager forcesManager(structure);

  // Get some example forces from D3
  Utils::Dftd3::Dftd3 d3;
  d3.initialize(structure, 1.0, 0.7875, 0.4289, 4.4407, Dftd3::Damping::BJ); // for PBE (does not actually matter)
  d3.calculate(Utils::Derivative::First);
  auto energy = d3.getEnergy();
  Utils::GradientCollection forces = -d3.getGradients();
  forces *= Utils::Constants::kCalPerMol_per_hartree; // conversion for more convenient numbers

  std::vector<Eigen::RowVector3d> originalForces;
  std::vector<Eigen::RowVector3d> originalInternalForces;

  // Check that the force transformation works
  for (int i = 0; i < forces.rows(); ++i) {
    originalForces.emplace_back(forces.row(i));
    auto internalForce = forcesManager.toInternalRepresentation(forces.row(i), i);
    originalInternalForces.push_back(internalForce);
    ASSERT_THAT(internalForce.norm(), DoubleNear(forces.row(i).norm(), 1e-6));
    auto globalForce = forcesManager.toGlobalRepresentation(internalForce, i);
    for (int j = 0; j < 3; ++j) {
      ASSERT_THAT(globalForce(j), DoubleNear(forces(i, j), 1e-6));
    }
  }

  // Rotate molecule (test 1)
  Eigen::RowVector3d translation = -Utils::Geometry::Properties::getCenterOfMass(structure);
  Utils::PositionCollection translatedPos =
      Utils::Geometry::Manipulations::translatePositions(structure.getPositions(), translation);
  Eigen::Vector3d rotationAxis(1.0, 1.0, 1.0);                      // axis is not important here
  Eigen::AngleAxisd angleAxis(M_PI / 3, rotationAxis.normalized()); // 60 degree rotation
  Eigen::Quaterniond rotation(angleAxis);
  Utils::PositionCollection rotatedPos(translatedPos.rows(), 3);
  for (int k = 0; k < translatedPos.rows(); k++) {
    rotatedPos.row(k) = rotation * translatedPos.row(k);
  }

  // Update forces manager
  forcesManager.modifyPositions(rotatedPos);

  // Calculate forces again:
  d3.initialize({structure.getElements(), rotatedPos}, 1.0, 0.7875, 0.4289, 4.4407, Dftd3::Damping::BJ);
  d3.calculate(Utils::Derivative::First);
  ASSERT_THAT(energy, DoubleNear(d3.getEnergy(), 1e-6)); // make sure energy is equal
  forces = -d3.getGradients() * Utils::Constants::kCalPerMol_per_hartree;

  for (int i = 0; i < forces.rows(); ++i) {
    auto internalForce = forcesManager.toInternalRepresentation(forces.row(i), i);
    for (int j = 0; j < 3; ++j) {
      ASSERT_FALSE(std::abs(forces(i, j) - originalForces.at(i)(j)) < 1e-6);
      ASSERT_THAT(internalForce(j), DoubleNear(originalInternalForces.at(i)(j), 1e-6));
    }
  }

  // Rotate molecule (test 2)
  Eigen::Vector3d rotationAxisNew(9.0, 3.0, 1.0);                       // different axis
  Eigen::AngleAxisd angleAxisNew(-M_PI / 4, rotationAxis.normalized()); // -45 degree rotation
  Eigen::Quaterniond rotationNew(angleAxisNew);
  for (int k = 0; k < translatedPos.rows(); k++) {
    rotatedPos.row(k) = rotationNew * translatedPos.row(k);
  }

  // Update forces manager
  forcesManager.modifyPositions(rotatedPos);

  // Calculate forces again:
  d3.initialize({structure.getElements(), rotatedPos}, 1.0, 0.7875, 0.4289, 4.4407, Dftd3::Damping::BJ);
  d3.calculate(Utils::Derivative::First);
  ASSERT_THAT(energy, DoubleNear(d3.getEnergy(), 1e-6)); // make sure energy is equal
  forces = -d3.getGradients() * Utils::Constants::kCalPerMol_per_hartree;

  for (int i = 0; i < forces.rows(); ++i) {
    auto internalForce = forcesManager.toInternalRepresentation(forces.row(i), i);
    for (int j = 0; j < 3; ++j) {
      ASSERT_FALSE(std::abs(forces(i, j) - originalForces.at(i)(j)) < 1e-6);
      ASSERT_THAT(internalForce(j), DoubleNear(originalInternalForces.at(i)(j), 1e-6));
    }
  }
}

TEST_F(AnAtomicForcesManagerTest, FeaturesAreCorrectlyGenerated) {
  std::stringstream ss1("13\n\n"
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
  auto structure = XyzStreamHandler::read(ss1);
  AtomicForcesManager forcesManager(structure);

  // For such small molecules, all atoms should be considered when creating the feature vector
  for (int i = 0; i < structure.size(); ++i) {
    auto features = forcesManager.calculateFeatures(i);
    ASSERT_THAT(features.size(), Eq(3 * (structure.size() - 1)));
  }

  // Feature vectors for atoms 0 and 1
  auto featuresZero = forcesManager.calculateFeatures(0);
  auto featuresOne = forcesManager.calculateFeatures(1);

  // Modify positions
  Utils::PositionCollection positions = structure.getPositions();
  // Push zero-th atom very far away
  for (int k = 0; k < 3; ++k) {
    positions(0, k) += 200.0;
  }
  forcesManager.modifyPositions(positions);

  // New feature vectors for atoms 0 and 1
  auto featuresZeroNew = forcesManager.calculateFeatures(0);
  auto featuresOneNew = forcesManager.calculateFeatures(1);

  // The size of the feature vectors should not change
  ASSERT_THAT(featuresZeroNew.size(), Eq(featuresZero.size()));
  ASSERT_THAT(featuresOneNew.size(), Eq(featuresOne.size()));

  for (int a = 0; a < featuresZero.size(); ++a) {
    ASSERT_TRUE(std::abs(featuresZero(a)) > std::abs(featuresZeroNew(a)));
    ASSERT_TRUE(std::abs(featuresZeroNew(a)) < 1e-4);
  }

  for (int b = 0; b < featuresOne.size(); ++b) {
    if ((featuresOneNew.size() - b) <= 3)
      ASSERT_TRUE(std::abs(featuresOneNew(b)) < 1e-4);
    else
      ASSERT_FALSE(std::abs(featuresOneNew(b)) < 1e-4);

    ASSERT_FALSE(std::abs(featuresOne(b)) < 1e-4);
  }

  // For large molecules, not all atoms should be considered when creating the feature vector
  std::stringstream ss2("148\n"
                        "\n"
                        "C   1.566858   -2.923359    1.426486 \n"
                        "C   1.934594   -2.154884    2.607139 \n"
                        "C   2.785736   -3.244685    0.697192 \n"
                        "C   0.403095   -2.607921    0.726967 \n"
                        "C   1.119268   -1.105634    3.038930 \n"
                        "C   3.384417   -2.000206    2.607052 \n"
                        "C   3.910629   -2.674581    1.424839 \n"
                        "C   2.785736   -3.244685   -0.697192 \n"
                        "C   0.403095   -2.607921   -0.726967 \n"
                        "C  -0.449292   -1.517912    1.177465 \n"
                        "C   1.719163    0.145401    3.487067 \n"
                        "C  -0.102196   -0.783264    2.311623 \n"
                        "C   3.957080   -0.801805    3.038093 \n"
                        "C   4.988594   -2.124083    0.728624 \n"
                        "C   1.566858   -2.923359   -1.426486 \n"
                        "C   3.910629   -2.674581   -1.424839 \n"
                        "C  -0.449292   -1.517912   -1.177465 \n"
                        "C  -0.981283   -0.844974    0.000000 \n"
                        "C   0.868417    1.241620    3.039715 \n"
                        "C   3.107735    0.293126    3.488645 \n"
                        "C  -0.257908    0.666843    2.312492 \n"
                        "C   5.081324   -0.227697    2.309235 \n"
                        "C   5.586099   -0.874428    1.180069 \n"
                        "C   4.988594   -2.124083   -0.728624 \n"
                        "C   1.934594   -2.154884   -2.607139 \n"
                        "C   3.384417   -2.000206   -2.607052 \n"
                        "C  -0.102196   -0.783264   -2.311623 \n"
                        "C  -1.133272    0.540745    0.000000 \n"
                        "C   1.442785    2.439068    2.606572 \n"
                        "C   3.707263    1.542787    3.038132 \n"
                        "C  -0.756403    1.309912    1.179654 \n"
                        "C   4.926738    1.221120    2.308793 \n"
                        "C   5.956118   -0.102109    0.000000 \n"
                        "C   5.586099   -0.874428   -1.180069 \n"
                        "C   1.119268   -1.105634   -3.038930 \n"
                        "C   3.957080   -0.801805   -3.038093 \n"
                        "C  -0.257908    0.666843   -2.312492 \n"
                        "C  -0.756403    1.309912   -1.179654 \n"
                        "C   2.892504    2.592196    2.607421 \n"
                        "C   0.920324    3.109489    1.422557 \n"
                        "C  -0.151213    2.553969    0.726292 \n"
                        "C   5.280777    1.961628    1.178786 \n"
                        "C   5.808001    1.285688    0.000000 \n"
                        "C   5.081324   -0.227697   -2.309235 \n"
                        "C   1.719163    0.145401   -3.487067 \n"
                        "C   3.107735    0.293126   -3.488645 \n"
                        "C   0.868417    1.241620   -3.039715 \n"
                        "C  -0.151213    2.553969   -0.726292 \n"
                        "C   3.262890    3.363236    1.427525 \n"
                        "C   2.044302    3.679589    0.696561 \n"
                        "C   4.432881    3.057001    0.729556 \n"
                        "C   5.280777    1.961628   -1.178786 \n"
                        "C   4.926738    1.221120   -2.308793 \n"
                        "C   3.707263    1.542787   -3.038132 \n"
                        "C   1.442785    2.439068   -2.606572 \n"
                        "C   0.920324    3.109489   -1.422557 \n"
                        "C   2.044302    3.679589   -0.696561 \n"
                        "C   4.432881    3.057001   -0.729556 \n"
                        "C   2.892504    2.592196   -2.607421 \n"
                        "C   3.262890    3.363236   -1.427525 \n"
                        "C  -4.134666    2.071846   -1.381300 \n"
                        "C  -3.404130    3.075861   -0.728345 \n"
                        "C  -3.404130    3.075861    0.728345 \n"
                        "C  -4.134666    2.071846    1.381300 \n"
                        "C  -4.869709    1.093881    0.706980 \n"
                        "C  -5.745992    0.202727    1.506195 \n"
                        "C  -6.683812    0.812836    2.357931 \n"
                        "C  -7.511269    0.059871    3.187984 \n"
                        "C  -7.399540   -1.331844    3.191066 \n"
                        "C  -6.463822   -1.948343    2.363837 \n"
                        "C  -5.634048   -1.201450    1.507502 \n"
                        "C  -4.626077   -1.940516    0.707658 \n"
                        "C  -3.736189   -2.779587    1.381906 \n"
                        "C  -2.840555   -3.639233    0.727981 \n"
                        "C  -2.840555   -3.639233   -0.727981 \n"
                        "C  -3.736189   -2.779587   -1.381906 \n"
                        "C  -4.626077   -1.940516   -0.707658 \n"
                        "C  -5.634048   -1.201450   -1.507502 \n"
                        "C  -6.463822   -1.948343   -2.363837 \n"
                        "C  -7.399540   -1.331844   -3.191066 \n"
                        "C  -7.511269    0.059871   -3.187984 \n"
                        "C  -6.683812    0.812836   -2.357931 \n"
                        "C  -5.745992    0.202727   -1.506195 \n"
                        "C  -4.869709    1.093881   -0.706980 \n"
                        "C  -2.610492    4.055256   -1.482051 \n"
                        "C  -2.124279    4.014626   -2.843111 \n"
                        "C  -1.056682    4.805336   -3.274392 \n"
                        "C  -0.365189    5.705042   -2.380159 \n"
                        "C   0.964976    6.275190   -2.450030 \n"
                        "C   1.626235    6.753961   -1.318866 \n"
                        "C   1.026857    6.714115    0.000000 \n"
                        "C   1.626235    6.753961    1.318866 \n"
                        "C   0.964976    6.275190    2.450030 \n"
                        "C  -0.365189    5.705042    2.380159 \n"
                        "C  -1.056682    4.805336    3.274392 \n"
                        "C  -2.124279    4.014626    2.843111 \n"
                        "C  -2.610492    4.055256    1.482051 \n"
                        "C  -2.061846    5.056332   -0.713522 \n"
                        "C  -0.977301    5.850840   -1.144956 \n"
                        "C  -0.304265    6.335072    0.000000 \n"
                        "C  -0.977301    5.850840    1.144956 \n"
                        "C  -2.061846    5.056332    0.713522 \n"
                        "C  -1.886630   -4.462192    1.480660 \n"
                        "C  -1.410215   -4.335713    2.839779 \n"
                        "C  -0.223636   -4.932798    3.271875 \n"
                        "C   0.608491   -5.708640    2.380659 \n"
                        "C   2.013561   -6.053094    2.450641 \n"
                        "C   2.745238   -6.415074    1.319296 \n"
                        "C   2.147627   -6.472791    0.000000 \n"
                        "C   2.745238   -6.415074   -1.319296 \n"
                        "C   2.013561   -6.053094   -2.450641 \n"
                        "C   0.608491   -5.708640   -2.380659 \n"
                        "C  -0.223636   -4.932798   -3.271875 \n"
                        "C  -1.410215   -4.335713   -2.839779 \n"
                        "C  -1.886630   -4.462192   -1.480660 \n"
                        "C  -1.174026   -5.356318    0.713596 \n"
                        "C   0.028569   -5.954463    1.145399 \n"
                        "C   0.771902   -6.321987    0.000000 \n"
                        "C   0.028569   -5.954463   -1.145399 \n"
                        "C  -1.174026   -5.356318   -0.713596 \n"
                        "H  -4.155476    2.066684   -2.474527 \n"
                        "H  -4.155476    2.066684    2.474527 \n"
                        "H  -6.760104    1.902454    2.350841 \n"
                        "H  -8.239479    0.557602    3.830662 \n"
                        "H  -8.039451   -1.936433    3.835300 \n"
                        "H  -6.367895   -3.036397    2.361827 \n"
                        "H  -3.757191   -2.778033    2.475074 \n"
                        "H  -3.757191   -2.778033   -2.475074 \n"
                        "H  -6.367895   -3.036397   -2.361827 \n"
                        "H  -8.039451   -1.936433   -3.835300 \n"
                        "H  -8.239479    0.557602   -3.830662 \n"
                        "H  -6.760104    1.902454   -2.350841 \n"
                        "H  -2.531759    3.279815   -3.541078 \n"
                        "H  -0.678519    4.658674   -4.289312 \n"
                        "H   1.512313    6.241491   -3.395986 \n"
                        "H   2.663909    7.081469   -1.424482 \n"
                        "H   2.663909    7.081469    1.424482 \n"
                        "H   1.512313    6.241491    3.395986 \n"
                        "H  -0.678519    4.658674    4.289312 \n"
                        "H  -2.531759    3.279815    3.541078 \n"
                        "H  -1.933288   -3.674655    3.535551 \n"
                        "H   0.123937   -4.721664    4.286617 \n"
                        "H   2.545261   -5.938178    3.399094 \n"
                        "H   3.822121   -6.569281    1.425436 \n"
                        "H   3.822121   -6.569281   -1.425436 \n"
                        "H   2.545261   -5.938178   -3.399094 \n"
                        "H   0.123937   -4.721664   -4.286617 \n"
                        "H  -1.933288   -3.674655   -3.535551 \n");
  auto largeSystem = XyzStreamHandler::read(ss2);
  AtomicForcesManager forcesManagerLarge(largeSystem);
  for (int index = 0; index < largeSystem.size(); ++index) {
    auto features = forcesManagerLarge.calculateFeatures(index);
    ASSERT_TRUE(features.size() < (3 * (largeSystem.size() - 1)));
  }
}

} // namespace Tests
} // namespace Scine
