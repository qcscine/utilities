/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/CalculatorBasics/CalculationRoutines.h>
#include <Utils/CalculatorBasics/TestCalculator.h>
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Utils {
namespace Tests {

class ACalculationRoutinesTest : public Test {
 public:
  TestCalculator calculator;
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

} // namespace Tests
} // namespace Utils
} // namespace Scine
