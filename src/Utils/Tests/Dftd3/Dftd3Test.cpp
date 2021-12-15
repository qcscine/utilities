/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Dftd3/Dftd3.h"
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/*
 * The usual namespace "Tests" is omitted, because the test classes
 * need to be declared friends of the Dftd3 class in Dftd3.h
 */
class ADftd3Calculation : public Test {
 public:
  ADftd3Calculation() {
  }
  Dftd3::Dftd3 method;
};

// Test, that the D3 energy is zero after the initialization if no calculation has been done yet.
TEST_F(ADftd3Calculation, EnergyIsZeroAfterInitialization) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "H     2.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 0.7875, 0.4289, 4.4407, Dftd3::Damping::BJ);
  EXPECT_THAT(method.getEnergy(), DoubleNear(0.0, 10e-10));
}

// Test, that the D3 gradient is zero after the initialization if no calculation has been done yet.
TEST_F(ADftd3Calculation, GradientIsZeroAfterInitialization) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "H     2.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 0.7875, 0.4289, 4.4407, Dftd3::Damping::BJ);
  auto gradient = method.getGradients();
  EXPECT_THAT(gradient(0), DoubleNear(0.0, 10e-8));
  EXPECT_THAT(gradient(1), DoubleNear(0.0, 10e-8));
  EXPECT_THAT(gradient(2), DoubleNear(0.0, 10e-8));
  EXPECT_THAT(gradient(3), DoubleNear(0.0, 10e-8));
  EXPECT_THAT(gradient(4), DoubleNear(0.0, 10e-8));
  EXPECT_THAT(gradient(5), DoubleNear(0.0, 10e-8));
}

// Test, that the correct number of atoms are stored after initialization.
TEST_F(ADftd3Calculation, HasCorrectNumberOfAtomsAfterInitialization) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "H     2.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 0.7875, 0.4289, 4.4407, Dftd3::Damping::BJ);
  ASSERT_THAT(static_cast<int>(method.getStructure().size()), Eq(2));
}

// Test, that the correct element types are present in the structure consisting of Dftd3::Dftd3Atoms.
TEST_F(ADftd3Calculation, HasCorrectAtomTypesAfterInitialization) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "H     2.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 0.7875, 0.4289, 4.4407, Dftd3::Damping::BJ);

  Dftd3::Dftd3Atom firstAtom = method.getStructure().at(0);
  Dftd3::Dftd3Atom secondAtom = method.getStructure().at(1);

  ASSERT_THAT(firstAtom.getElementType(), Eq(ElementType::C));
  ASSERT_THAT(secondAtom.getElementType(), Eq(ElementType::H));
}

// Test, that the atom indices are stored correctly in the structure consisting of Dftd3::Dftd3Atoms.
TEST_F(ADftd3Calculation, HasCorrectAtomIndecesAfterInitialization) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "H     2.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 0.7875, 0.4289, 4.4407, Dftd3::Damping::BJ);

  Dftd3::Dftd3Atom firstAtom = method.getStructure().at(0);
  Dftd3::Dftd3Atom secondAtom = method.getStructure().at(1);

  ASSERT_THAT(firstAtom.getIndex(), Eq(0));
  ASSERT_THAT(secondAtom.getIndex(), Eq(1));
}

// Test, that the correct coordination numbers are calculated.
TEST_F(ADftd3Calculation, HasCorrectCoordinationNumbers) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "H     2.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 0.7875, 0.4289, 4.4407, Dftd3::Damping::BJ);
  method.calculate(Derivative::None);

  Dftd3::Dftd3Atom firstAtom = method.getStructure().at(0);
  Dftd3::Dftd3Atom secondAtom = method.getStructure().at(1);

  // Correct coordination numbers were obtained with Grimme's implementation of D3
  EXPECT_THAT(firstAtom.getCoordinationNumber(), DoubleNear(0.010, 10e-4));
  EXPECT_THAT(secondAtom.getCoordinationNumber(), DoubleNear(0.010, 10e-4));
}

// Test, that the correct energy and gradients are calculated for the CO molecule (PBE functional).
TEST_F(ADftd3Calculation, HasCorrectDispersionEnergyandGradientForPBEWithBJDamping) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "O     1.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 0.7875, 0.4289, 4.4407, Dftd3::Damping::BJ); // Parameters for PBE
  method.calculate(Derivative::First);
  auto gradient = method.getGradients();

  // Correct energy and gradient were obtained with Grimme's implementation of D3
  EXPECT_THAT(method.getEnergy(), DoubleNear(-0.00041108, 10e-8));
  EXPECT_THAT(gradient.rows(), Eq(2));
  EXPECT_THAT(gradient.size(), Eq(6));
  EXPECT_THAT(gradient(0), DoubleNear(-0.5466834209e-6, 10e-8));
  EXPECT_THAT(gradient(1), DoubleNear(0.0, 10e-8));
  EXPECT_THAT(gradient(2), DoubleNear(0.0, 10e-8));
  EXPECT_THAT(gradient(3), DoubleNear(0.5466834209e-6, 10e-8));
  EXPECT_THAT(gradient(4), DoubleNear(0.0, 10e-8));
  EXPECT_THAT(gradient(5), DoubleNear(0.0, 10e-8));
}

TEST_F(ADftd3Calculation, HasCorrectDispersionEnergyandGradientForPBEWithZeroDamping) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "O     1.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 0.7220, 1.217, 14, Dftd3::Damping::Zero); // Parameters for PBE
  method.calculate(Derivative::First);
  auto gradient = method.getGradients();

  // Correct energy and gradient were obtained with Grimme's implementation of D3
  EXPECT_THAT(method.getEnergy(), DoubleNear(-0.00000012, 10e-8));
  EXPECT_THAT(gradient.rows(), Eq(2));
  EXPECT_THAT(gradient.size(), Eq(6));
  EXPECT_THAT(gradient(0), DoubleNear(0.51108861161807e-6, 10e-8));
  EXPECT_THAT(gradient(1), DoubleNear(0.0, 10e-8));
  EXPECT_THAT(gradient(2), DoubleNear(0.0, 10e-8));
  EXPECT_THAT(gradient(3), DoubleNear(-0.51108861161807e-6, 10e-8));
  EXPECT_THAT(gradient(4), DoubleNear(0.0, 10e-8));
  EXPECT_THAT(gradient(5), DoubleNear(0.0, 10e-8));
}

// Test, that the correct energy is calculated for the CO molecule (BP86 functional).
TEST_F(ADftd3Calculation, HasCorrectDispersionEnergyForBP86WithBJDamping) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "O     1.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 3.2822, 0.3946, 4.8516, Dftd3::Damping::BJ); // Parameters for BP86
  method.calculate(Derivative::None);

  // Correct energy was obtained with Grimme's implementation of D3
  EXPECT_THAT(method.getEnergy(), DoubleNear(-0.00062281, 10e-8));
}

TEST_F(ADftd3Calculation, HasCorrectDispersionEnergyForBP86WithZeroDamping) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "O     1.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 1.6830, 1.1390, 14, Dftd3::Damping::Zero); // Parameters for BP86
  method.calculate(Derivative::None);

  // Correct energy was obtained with Grimme's implementation of D3
  EXPECT_THAT(method.getEnergy(), DoubleNear(-0.00000028, 10e-8));
}

// Test, that the correct energy is calculated for the CO2 molecule (TPSS functional).
TEST_F(ADftd3Calculation, HasCorrectDispersionEnergyandGradientForTPSSWithBJDamping) {
  std::stringstream ss("3\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "O     1.0000000000    0.0000000000    0.0000000000\n"
                       "O     -1.000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 1.9435, 0.4535, 4.4752, Dftd3::Damping::BJ); // Parameters for TPSS
  method.calculate(Derivative::First);

  // Correct values were obtained with Grimme's implementation of D3
  EXPECT_THAT(method.getEnergy(), DoubleNear(-0.00116406, 10e-8));
  EXPECT_THAT(method.getGradients()(0), DoubleNear(0.0, 10e-8));
  EXPECT_THAT(method.getGradients()(3), DoubleNear(0.91470925678222e-6, 10e-8));
  EXPECT_THAT(method.getGradients()(6), DoubleNear(-0.91470925678222e-6, 10e-8));
}

TEST_F(ADftd3Calculation, HasCorrectDispersionEnergyandGradientForTPSSWithZeroDamping) {
  std::stringstream ss("3\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "O     1.0000000000    0.0000000000    0.0000000000\n"
                       "O     -1.000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 1.1050, 1.1660, 14, Dftd3::Damping::Zero); // Parameters for TPSS
  method.calculate(Derivative::First);

  // Correct values were obtained with Grimme's implementation of D3
  EXPECT_THAT(method.getEnergy(), DoubleNear(-0.00004030, 10e-8));
  EXPECT_THAT(method.getGradients()(0), DoubleNear(0.0, 10e-8));
  EXPECT_THAT(method.getGradients()(3), DoubleNear(-0.85163670605709e-4, 10e-8));
  EXPECT_THAT(method.getGradients()(6), DoubleNear(0.85163670605709e-4, 10e-8));
}

// Test, that the total energy calculated for a molecule with two atoms is equal to the energy calculated
// for just this one atom pair.
TEST_F(ADftd3Calculation, EnergyCalculationIsInternallyConsistent) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "O     1.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 1.9435, 0.4535, 4.4752, Dftd3::Damping::BJ); // Parameters for TPSS
  method.calculate(Derivative::None);
  double energy1 = method.getEnergy();

  Dftd3::Dftd3Atom firstAtom = method.getStructure().at(0);
  Dftd3::Dftd3Atom secondAtom = method.getStructure().at(1);
  double energy2 = method.evaluateEnergy(firstAtom, secondAtom);

  EXPECT_THAT(energy1, DoubleNear(energy2, 10e-8));
}

// Test, that the correct C6 coefficient is calculated for the CO molecule.
TEST_F(ADftd3Calculation, HasCorrectC6Coefficients) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "O     1.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 1.9435, 0.4535, 4.4752, Dftd3::Damping::BJ); // Parameters for TPSS
  method.calculate(Derivative::None);

  Dftd3::Dftd3Atom firstAtom = method.getStructure().at(0);
  Dftd3::Dftd3Atom secondAtom = method.getStructure().at(1);
  double c6 = method.calculateC6Coefficient(firstAtom, secondAtom);

  // Correct C6 coefficient were obtained with Grimme's implementation of D3
  EXPECT_THAT(c6, DoubleNear(0.224714e2, 10e-5));
}

// Test, that the derivatives of the C6 coefficients w.r.t the coordination numbers are correctly calculated.
TEST_F(ADftd3Calculation, HasCorrectDerivativeOfC6WrtCoordNumber) {
  std::stringstream ss("3\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "O     1.0000000000    0.0000000000    0.0000000000\n"
                       "O     -1.000000000    0.0000000000    0.0000000000\n");
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 1.9435, 0.4535, 4.4752, Dftd3::Damping::BJ); // Parameters for TPSS
  Dftd3::Dftd3Atom firstAtom = method.getStructure().at(0);
  Dftd3::Dftd3Atom secondAtom = method.getStructure().at(1);
  // Set the coordination numbers.
  firstAtom.setCoordinationNumber(1.5);
  secondAtom.setCoordinationNumber(2.775);
  // Do a numerical derivative.
  auto c6Old = method.calculateC6Coefficient(firstAtom, secondAtom);
  auto grad = method.evaluateGradientOfC6WrtCoordNumber(firstAtom, secondAtom);
  double stepsize = 0.00001;
  firstAtom.setCoordinationNumber(1.5 + stepsize);
  auto c6New = method.calculateC6Coefficient(firstAtom, secondAtom);
  auto numGrad = (c6New - c6Old) / stepsize;
  // Check the gradient calculation against the numerical gradient.
  EXPECT_THAT(grad, DoubleNear(numGrad, 10e-6));
}

// Test, that the partial derivative of the energy w.r.t to the C-O distance is correctly calculated.
// The test is performed against a numerical derivative.
TEST_F(ADftd3Calculation, HasCorrectDerivativeOfEnergyWrtDistance) {
  std::stringstream ss1("2\n\n"
                        "C     0.0000000000    0.0000000000   -0.0000000000\n"
                        "O     1.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection1 = XyzStreamHandler::read(ss1);
  method.initialize(atomCollection1, 1.0, 1.9435, 0.4535, 4.4752, Dftd3::Damping::BJ);

  method.calculateValuesForConstants();
  auto firstAtom = method.getStructure()[0];
  auto secondAtom = method.getStructure()[1];

  auto grad = method.evaluateGradientsWrtDistances(firstAtom, secondAtom);
  auto energy1 = method.evaluateEnergy(firstAtom, secondAtom);

  std::stringstream ss2("2\n\n"
                        "C     0.0000100000    0.0000000000   -0.0000000000\n"
                        "O     1.0000000000    0.0000000000    0.0000000000\n");
  auto atomCollection2 = XyzStreamHandler::read(ss2);

  method.initialize(atomCollection2, 1.0, 1.9435, 0.4535, 4.4752, Dftd3::Damping::BJ);
  method.calculateValuesForConstants();
  auto firstAtomNew = method.getStructure()[0];
  auto secondAtomNew = method.getStructure()[1];

  auto energy2 = method.evaluateEnergy(firstAtomNew, secondAtomNew);

  auto numGrad = (energy1 - energy2) / (0.00001 / (Constants::bohrRadius * Constants::angstrom_per_meter));

  EXPECT_THAT(grad, DoubleNear(numGrad, 10e-11));
}

// Test, that the correct D3 energy and gradients are calculated for a large molecule (PBE0 functional).
TEST_F(ADftd3Calculation, HasCorrectDispersionEnergyandGradientForPBE0WithBJDamping) {
  std::stringstream ss("148\n"
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
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 1.2177, 0.4145, 4.8593, Dftd3::Damping::BJ); // Parameters for PBE0
  method.calculate(Derivative::First);
  auto gradient = method.getGradients();

  // Correct energy and gradient were obtained with Grimme's implementation of D3
  EXPECT_THAT(method.getEnergy(), DoubleNear(-0.39503998, 10e-8));
  EXPECT_THAT(gradient(0), DoubleNear(0.13310894565617e-2, 10e-8));
  EXPECT_THAT(gradient(1), DoubleNear(0.12174773980560e-2, 10e-8));
  EXPECT_THAT(gradient(2), DoubleNear(0.12324240976931e-3, 10e-8));
  EXPECT_THAT(gradient(9), DoubleNear(0.44850151276792e-4, 10e-8));
  EXPECT_THAT(gradient(10), DoubleNear(0.18403297368559e-2, 10e-8));
  EXPECT_THAT(gradient(18), DoubleNear(-0.81163403700319e-3, 10e-8));
}

TEST_F(ADftd3Calculation, HasCorrectDispersionEnergyandGradientForPBE0WithZeroDamping) {
  std::stringstream ss("148\n"
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
  auto atomCollection = XyzStreamHandler::read(ss);
  method.initialize(atomCollection, 1.0, 0.9280, 1.2870, 14.0, Dftd3::Damping::Zero); // Parameters for PBE0
  method.calculate(Derivative::First);
  auto gradient = method.getGradients();

  // Correct energy and gradient were obtained with Grimme's implementation of D3
  EXPECT_THAT(method.getEnergy(), DoubleNear(-0.20543164, 10e-8));
  EXPECT_THAT(gradient(0), DoubleNear(0.62507533604589e-3, 10e-8));
  EXPECT_THAT(gradient(1), DoubleNear(0.68043605384560e-3, 10e-8));
  EXPECT_THAT(gradient(2), DoubleNear(0.67020089194758e-4, 10e-8));
  EXPECT_THAT(gradient(9), DoubleNear(0.64052000470617e-4, 10e-8));
  EXPECT_THAT(gradient(10), DoubleNear(0.97337076868985e-3, 10e-8));
  EXPECT_THAT(gradient(18), DoubleNear(-0.26099305115032e-3, 10e-8));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
