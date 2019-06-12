/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Constants.h>
#include <Utils/Geometry/ElementInfo.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class ElementInfoTest ElementInfoTest.cpp
 * @brief Comprises tests for the class Scine::Utils::ElementInfo.
 * @test
 */
class ElementInfoTest : public Test {
 public:
 protected:
  void SetUp() override {
  }
};

TEST_F(ElementInfoTest, CanConvertSymbolStringToElementType) {
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("H"), Eq(ElementType::H));
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("C"), Eq(ElementType::C));
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("Rn"), Eq(ElementType::Rn));
}

TEST_F(ElementInfoTest, CanConvertElementTypeToSymbol) {
  ASSERT_THAT(ElementInfo::symbol(ElementType::H), Eq(std::string("H")));
  ASSERT_THAT(ElementInfo::symbol(ElementType::C), Eq(std::string("C")));
  ASSERT_THAT(ElementInfo::symbol(ElementType::Rn), Eq(std::string("Rn")));
}

TEST_F(ElementInfoTest, CanGetVdwRadiusForElementInAtomicUnits) {
  double vdwInPicometers = 210.0;
  double vdwInAngstrom = vdwInPicometers / 100.0;
  double vdwInAtomicUnits = vdwInAngstrom * Constants::bohr_per_angstrom;
  ASSERT_THAT(ElementInfo::vdwRadius(ElementType::Si), Eq(vdwInAtomicUnits));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
