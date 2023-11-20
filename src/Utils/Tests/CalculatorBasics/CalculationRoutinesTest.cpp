/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/CalculatorBasics.h>
#include <Utils/Geometry/AtomCollection.h>
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Utils {
namespace Tests {

class ACalculationRoutinesTest : public Test {
 public:
  TestCalculator calculator;
  ElementTypeCollection arbitraryElements;
  PositionCollection randomPositions;
  AtomCollection randomAtomCollection;

 private:
  void SetUp() final {
    arbitraryElements = ElementTypeCollection{ElementType::Am, ElementType::Ce, ElementType::Ca};
    randomPositions = Eigen::MatrixX3d::Random(3, 3);
    randomAtomCollection.resize(arbitraryElements.size());
    randomAtomCollection.setElements(arbitraryElements);
    randomAtomCollection.setPositions(randomPositions);
  }
};

TEST_F(ACalculationRoutinesTest, CanSetLog) {
  CalculationRoutines::setLog(calculator, false, false, false);
  auto log = calculator.getLog();
  ASSERT_FALSE(log.error);
  ASSERT_FALSE(log.warning);
  ASSERT_FALSE(log.output);

  CalculationRoutines::setLog(calculator, true, true, true);
  log = calculator.getLog();
  ASSERT_TRUE(log.error);
  ASSERT_TRUE(log.warning);
  ASSERT_TRUE(log.output);

  CalculationRoutines::setLog(calculator, true, false, true);
  log = calculator.getLog();
  ASSERT_TRUE(log.error);
  ASSERT_FALSE(log.warning);
  ASSERT_TRUE(log.output);

  CalculationRoutines::setLog(calculator, true, false, false);
  log = calculator.getLog();
  ASSERT_TRUE(log.error);
  ASSERT_FALSE(log.warning);
  ASSERT_FALSE(log.output);
}

TEST_F(ACalculationRoutinesTest, SpinPropensity) {
  calculator.setStructure(randomAtomCollection);
  auto refSettings = calculator.settings();
  auto log = Core::Log::silent();
  // 0 does nothing
  auto calc = CalculationRoutines::spinPropensity(calculator, log, 0);
  ASSERT_THAT(calc->settings().getInt(SettingsNames::spinMultiplicity), 1);
  // actual check
  int range = 2;
  calc = CalculationRoutines::spinPropensity(calculator, log, range);
  ASSERT_THAT(calc->settings().getInt(SettingsNames::spinMultiplicity), 2 * range + 1);
}

TEST_F(ACalculationRoutinesTest, DispersionSplitting) {
  std::string correct = "pbe-d3bj";
  auto outCorrect = CalculationRoutines::splitIntoMethodAndDispersion(correct);
  ASSERT_THAT(outCorrect.first, std::string{"pbe"});
  ASSERT_THAT(outCorrect.second, std::string{"d3bj"});

  std::string correctException = "dlpno-ccsd(t)";
  auto outCorrectException = CalculationRoutines::splitIntoMethodAndDispersion(correctException);
  ASSERT_THAT(outCorrectException.first, std::string{"dlpno-ccsd(t)"});
  ASSERT_THAT(outCorrectException.second, std::string{""});

  std::string noDisp = "pbe";
  auto outNoDisp = CalculationRoutines::splitIntoMethodAndDispersion(noDisp);
  ASSERT_THAT(outNoDisp.first, std::string{"pbe"});
  ASSERT_TRUE(outNoDisp.second.empty());

  std::string tooMuch = "pbe-d3-bj";
  ASSERT_THROW(CalculationRoutines::splitIntoMethodAndDispersion(tooMuch), std::logic_error);

  std::string hackyInput = "pbe special_scf_setting";
  ASSERT_THROW(CalculationRoutines::splitIntoMethodAndDispersion(hackyInput), std::logic_error);
  hackyInput += "-d3bj";
  ASSERT_THROW(CalculationRoutines::splitIntoMethodAndDispersion(hackyInput), std::logic_error);
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
