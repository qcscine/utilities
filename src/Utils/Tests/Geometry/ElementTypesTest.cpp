/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry/ElementTypes.h"
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class ElementTypeTest ElementTypesTest.cpp
 * @brief Comprises tests for the class Scine::Utils::ElementTypes.
 * @test
 */
class ElementTypeTest : public Test {
 public:
  ElementType H_type = ElementType::H;
  ElementType C_type = ElementType::C;
  ElementType Rn_type = ElementType::Rn;
  int H_Z = 1;
  int C_Z = 6;
  int Rn_Z = 86;

 protected:
  void SetUp() override {
  }
};

TEST_F(ElementTypeTest, HasImplicitConversionFromAtomicNumber) {
  ASSERT_THAT(static_cast<ElementType>(H_Z), Eq(H_type));
  ASSERT_THAT(static_cast<ElementType>(C_Z), Eq(C_type));
  ASSERT_THAT(static_cast<ElementType>(Rn_Z), Eq(Rn_type));
}

TEST_F(ElementTypeTest, HasImplicitConversionToAtomicNumber) {
  ASSERT_THAT(static_cast<int>(H_type), Eq(H_Z));
  ASSERT_THAT(static_cast<int>(C_type), Eq(C_Z));
  ASSERT_THAT(static_cast<int>(Rn_type), Eq(Rn_Z));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
