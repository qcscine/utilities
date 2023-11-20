/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Scf/LcaoUtils/SpinMode.h"
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class Scine::Utils::Tests::SpinModeTests
 * @brief Comprises tests for the class Scine::Utils::SpinModeInterpreter.
 * @test
 */
TEST(SpinModeTests, isSelfConsistent) {
  SpinMode restricted = SpinMode::Restricted;
  SpinMode unrestricted = SpinMode::Unrestricted;
  SpinMode restrictedOpenShell = SpinMode::RestrictedOpenShell;
  SpinMode any = SpinMode::Any;
  SpinMode none = SpinMode::None;

  EXPECT_EQ("restricted", SpinModeInterpreter::getStringFromSpinMode(restricted));
  EXPECT_EQ("unrestricted", SpinModeInterpreter::getStringFromSpinMode(unrestricted));
  EXPECT_EQ("restricted_open_shell", SpinModeInterpreter::getStringFromSpinMode(restrictedOpenShell));
  EXPECT_EQ("any", SpinModeInterpreter::getStringFromSpinMode(any));
  EXPECT_EQ("none", SpinModeInterpreter::getStringFromSpinMode(none));

  EXPECT_EQ(restricted, SpinModeInterpreter::getSpinModeFromString("restricted"));
  EXPECT_EQ(unrestricted, SpinModeInterpreter::getSpinModeFromString("unrestricted"));
  EXPECT_EQ(restrictedOpenShell, SpinModeInterpreter::getSpinModeFromString("restricted_open_shell"));
  EXPECT_EQ(any, SpinModeInterpreter::getSpinModeFromString("any"));
  EXPECT_EQ(none, SpinModeInterpreter::getSpinModeFromString("none"));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
