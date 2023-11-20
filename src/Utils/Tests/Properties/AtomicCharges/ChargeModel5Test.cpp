/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Properties/AtomicCharges/ChargeModel5.h>
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Utils {

/*
 * The tests are based on the results from the paper J. Chem. Theory Comput. 2012, 8, 2, 527-541
 * found in the Supporting Information.
 */

/**
 * @class AChargeModel5Test ChargeModel5Test.cpp
 * @brief Tests the Charge Model 5 corrections to the Hirshfeld population analysis.
 * @test
 */
class AChargeModel5Test : public Test {};

TEST_F(AChargeModel5Test, Cm5ChargesAreCorrectlyCalculatedForWater) {
  std::stringstream ss("3\n\n"
                       "O       0.000000    0.000000    0.117551\n"
                       "H       0.000000    0.756673   -0.470206\n"
                       "H       0.000000   -0.756673   -0.470206\n");
  auto structure = Utils::XyzStreamHandler::read(ss);

  std::vector<double> hirshfeldChargesWrongSize = {-0.3065, 0.1533};
  EXPECT_THROW(ChargeModel5::calculateCm5Charges(hirshfeldChargesWrongSize, structure), std::runtime_error);

  std::vector<double> hirshfeldCharges = {-0.3065, 0.1533, 0.1533};
  std::vector<double> cm5Charges = ChargeModel5::calculateCm5Charges(hirshfeldCharges, structure);

  ASSERT_THAT(cm5Charges.size(), Eq(3));
  ASSERT_THAT(cm5Charges.at(0), DoubleNear(-0.6423, 1e-4));
  ASSERT_THAT(cm5Charges.at(1), DoubleNear(0.3211, 1e-4));
  ASSERT_THAT(cm5Charges.at(2), DoubleNear(0.3211, 1e-4));
}

TEST_F(AChargeModel5Test, Cm5ChargesAreCorrectlyCalculatedFor18AtomMolecule) {
  std::stringstream ss("18\n\n"
                       "C      -2.585060   -1.031121   -0.039836\n"
                       "C      -3.019151    0.286906    0.001659\n"
                       "C      -2.077282    1.304388    0.048828\n"
                       "C      -0.716691    1.033287    0.054999\n"
                       "C      -0.285145   -0.293334    0.012905\n"
                       "C      -1.232035   -1.319904   -0.034006\n"
                       "N       1.063351   -0.670447    0.030124\n"
                       "C       2.181849    0.137092   -0.019814\n"
                       "N       3.362904   -0.587503   -0.014599\n"
                       "O       2.161817    1.345506   -0.105228\n"
                       "H       1.221439   -1.662067   -0.028326\n"
                       "H       4.173941   -0.013589    0.137055\n"
                       "H       3.393655   -1.474935    0.457432\n"
                       "H       0.009797    1.828633    0.085145\n"
                       "H      -0.897370   -2.351338   -0.066181\n"
                       "H      -2.401084    2.337170    0.080642\n"
                       "H      -4.076370    0.515908   -0.003521\n"
                       "H      -3.301242   -1.842064   -0.077515\n");

  auto structure = Utils::XyzStreamHandler::read(ss);

  std::vector<double> hirshfeldCharges = {-0.0393, -0.0511, -0.0335, -0.0533, 0.0522, -0.0630, -0.0866, 0.2076, -0.1753,
                                          -0.3077, 0.1241,  0.1299,  0.1187,  0.0340, 0.0304,  0.0394,  0.0356, 0.0379};
  std::vector<double> correctCm5Charges = {-0.0950, -0.1071, -0.0891, -0.1016, 0.1053, -0.1113,
                                           -0.4302, 0.3639,  -0.6233, -0.3517, 0.3272, 0.3259,
                                           0.3139,  0.1017,  0.0900,  0.0956,  0.0917, 0.0941};
  std::vector<double> cm5Charges = ChargeModel5::calculateCm5Charges(hirshfeldCharges, structure);

  ASSERT_THAT(cm5Charges.size(), Eq(structure.size()));
  for (int i = 0; i < structure.size(); ++i) {
    ASSERT_THAT(cm5Charges.at(i), DoubleNear(correctCm5Charges.at(i), 1e-4));
  }
}

TEST_F(AChargeModel5Test, Cm5ChargesAreCorrectlyCalculatedForUnusualMolecule) {
  std::stringstream ss("6\n\n"
                       "C      -0.609205   -0.729240    0.000221\n"
                       "C      -1.930996   -0.681585    0.000231\n"
                       "Na      1.673943   -0.877925    0.000017\n"
                       "Cl      0.107604    1.151032   -0.000156\n"
                       "H      -2.564230    0.197873   -0.000122\n"
                       "H      -2.437197   -1.643297   -0.000112\n");

  auto structure = Utils::XyzStreamHandler::read(ss);

  std::vector<double> hirshfeldCharges = {-0.2004, -0.1366, 0.5685, -0.2929, 0.0314, 0.0299};
  std::vector<double> correctCm5Charges = {-0.2818, -0.2361, 0.7136, -0.3628, 0.0845, 0.0827};
  std::vector<double> cm5Charges = ChargeModel5::calculateCm5Charges(hirshfeldCharges, structure);

  ASSERT_THAT(cm5Charges.size(), Eq(structure.size()));
  for (int i = 0; i < structure.size(); ++i) {
    ASSERT_THAT(cm5Charges.at(i), DoubleNear(correctCm5Charges.at(i), 1e-4));
  }
}

TEST_F(AChargeModel5Test, Cm5ChargesAreCorrectlyCalculatedFor16AtomMolecule) {
  std::stringstream ss("16\n\n"
                       "C       0.000000    0.762632   -1.409539\n"
                       "C       0.000000   -0.762632   -1.409539\n"
                       "C      -0.518994   -1.107609   -0.015440\n"
                       "C       0.000000    0.000000    0.870373\n"
                       "C       0.518994    1.107609   -0.015440\n"
                       "H       0.591966    1.199648   -2.212395\n"
                       "H      -1.022171    1.130592   -1.529771\n"
                       "H       1.022171   -1.130592   -1.529771\n"
                       "H      -0.591966   -1.199648   -2.212395\n"
                       "H      -1.613303   -1.086968   -0.014087\n"
                       "H      -0.231847   -2.098777    0.334583\n"
                       "H       0.231847    2.098777    0.334583\n"
                       "H       1.613303    1.086968   -0.014087\n"
                       "C       0.000000    0.000000    2.197574\n"
                       "H       0.372990    0.841577    2.767704\n"
                       "H      -0.372990   -0.841577    2.767704\n");

  auto structure = Utils::XyzStreamHandler::read(ss);

  std::vector<double> hirshfeldCharges = {-0.0484, -0.0484, -0.0483, 0.0101, -0.0483, 0.0286,  0.0264, 0.0264,
                                          0.0286,  0.0328,  0.0333,  0.0333, 0.0328,  -0.1148, 0.0279, 0.0279};
  std::vector<double> correctCm5Charges = {-0.1575, -0.1575, -0.1521, -0.0128, -0.1521, 0.0829,  0.0821, 0.0821,
                                           0.0829,  0.0886,  0.0884,  0.0884,  0.0886,  -0.2145, 0.0812, 0.0812};
  std::vector<double> cm5Charges = ChargeModel5::calculateCm5Charges(hirshfeldCharges, structure);

  ASSERT_THAT(cm5Charges.size(), Eq(structure.size()));
  for (int i = 0; i < structure.size(); ++i) {
    ASSERT_THAT(cm5Charges.at(i), DoubleNear(correctCm5Charges.at(i), 1e-4));
  }
}

TEST_F(AChargeModel5Test, Cm5ChargesAreCorrectlyCalculatedForHydrogenChloride) {
  std::stringstream ss("2\n\n"
                       "Cl      0.000000    0.000000    0.070785\n"
                       "H       0.000000    0.000000   -1.203341\n");

  auto structure = Utils::XyzStreamHandler::read(ss);

  std::vector<double> hirshfeldCharges = {-0.1358, 0.1358};
  std::vector<double> correctCm5Charges = {-0.1918, 0.1918};
  std::vector<double> cm5Charges = ChargeModel5::calculateCm5Charges(hirshfeldCharges, structure);

  ASSERT_THAT(cm5Charges.size(), Eq(structure.size()));
  for (int i = 0; i < structure.size(); ++i) {
    ASSERT_THAT(cm5Charges.at(i), DoubleNear(correctCm5Charges.at(i), 1e-4));
  }
}

} // namespace Utils
} // namespace Scine
