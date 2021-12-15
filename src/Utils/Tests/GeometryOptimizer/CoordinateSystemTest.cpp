/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/GeometryOptimization/CoordinateSystem.h"
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class Scine::Utils::Tests::CoordinateSystemTests
 * @brief Comprises tests for the class Scine::Utils::CoordinateSystemInterpreter.
 * @test
 */
TEST(CoordinateSystemTests, isSelfConsistent) {
  CoordinateSystem cartesian = CoordinateSystem::Cartesian;
  CoordinateSystem cartesianWithoutRotTrans = CoordinateSystem::CartesianWithoutRotTrans;
  CoordinateSystem internal = CoordinateSystem::Internal;

  EXPECT_EQ("cartesian", CoordinateSystemInterpreter::getStringFromCoordinateSystem(cartesian));
  EXPECT_EQ("cartesianWithoutRotTrans", CoordinateSystemInterpreter::getStringFromCoordinateSystem(cartesianWithoutRotTrans));
  EXPECT_EQ("internal", CoordinateSystemInterpreter::getStringFromCoordinateSystem(internal));

  EXPECT_EQ(cartesian, CoordinateSystemInterpreter::getCoordinateSystemFromString("cartesian"));
  EXPECT_EQ(cartesianWithoutRotTrans,
            CoordinateSystemInterpreter::getCoordinateSystemFromString("cartesianWithoutRotTrans"));
  EXPECT_EQ(internal, CoordinateSystemInterpreter::getCoordinateSystemFromString("internal"));

  EXPECT_EQ(cartesian, CoordinateSystemInterpreter::getCoordinateSystemFromString(
                           CoordinateSystemInterpreter::getStringFromCoordinateSystem(cartesian)));
  EXPECT_EQ(cartesianWithoutRotTrans,
            CoordinateSystemInterpreter::getCoordinateSystemFromString(
                CoordinateSystemInterpreter::getStringFromCoordinateSystem(cartesianWithoutRotTrans)));
  EXPECT_EQ(internal, CoordinateSystemInterpreter::getCoordinateSystemFromString(
                          CoordinateSystemInterpreter::getStringFromCoordinateSystem(internal)));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
