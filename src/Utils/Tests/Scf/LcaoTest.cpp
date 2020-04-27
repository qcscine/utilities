/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Scf/MethodInterfaces/LcaoMethod.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

class LcaoTestMethod : public LcaoMethod {
 public:
  LcaoTestMethod() : LcaoMethod(true, Utils::derivOrder::two, true) {
  }
  void initialize() final {
    molecularCharge_ = 0;
    nElectronsForUnchargedSpecies_ = 10;
    nAOs_ = 10;
  }
};

/**
 * @class ALcaoTest LcaoTest.cpp
 * @brief Comprises tests for the class Scine::Utils::LcaoMethod.
 * @test
 */
class ALcaoTest : public Test {
 public:
  LcaoTestMethod lcaoMethod_;

 protected:
  void SetUp() override {
    lcaoMethod_.initialize();
  }
};

TEST_F(ALcaoTest, CorrectlyChangesUnrestrictedCalculation) {
  lcaoMethod_.setMolecularCharge(-2);
  lcaoMethod_.setSpinMultiplicity(3);
  lcaoMethod_.setUnrestrictedCalculation(false);

  lcaoMethod_.verifyPesValidity();

  ASSERT_THAT(lcaoMethod_.getMolecularCharge(), -2);
  ASSERT_THAT(lcaoMethod_.unrestrictedCalculationRunning(), true);
  ASSERT_THAT(lcaoMethod_.spinMultiplicity(), 3);
}

TEST_F(ALcaoTest, CorrectlyChangesSpinMultiplicity) {
  lcaoMethod_.setMolecularCharge(-2);
  lcaoMethod_.setSpinMultiplicity(2);
  lcaoMethod_.setUnrestrictedCalculation(true);

  lcaoMethod_.verifyPesValidity();

  ASSERT_THAT(lcaoMethod_.getMolecularCharge(), -2);
  ASSERT_THAT(lcaoMethod_.unrestrictedCalculationRunning(), true);
  ASSERT_THAT(lcaoMethod_.spinMultiplicity(), 1);
}

// TODO: Add more unit tests

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
