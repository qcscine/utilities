/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/BasisSet.h>
#include <Utils/DataStructures/IntegralSpecifier.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Utils {
namespace Integrals {

using namespace testing;

class IntegralSpecifierTest : public Test {};
class BasisSetTest : public Test {};

TEST_F(IntegralSpecifierTest, CanGetCorrectParticleType) {
  auto type = getParticleType("e");
  ASSERT_EQ(type.charge, -1);
  ASSERT_EQ(type.mass, 1);
  ASSERT_EQ(type.spin, 1);
  ASSERT_EQ(type.symbol, Utils::ElementType::E);
}

TEST_F(BasisSetTest, CanCorrectlyBeConstructed) {
  AtomCollection atoms;
  std::string name = "test";
  BasisSet basis(name, atoms);
}

TEST_F(BasisSetTest, CorrectIDs) {
  AtomCollection atoms;
  std::string name = "test";
  BasisSet basis(name, atoms);

  auto basis2 = basis;

  auto basis3(basis);

  ASSERT_EQ(basis.getID(), basis2.getID());

  ASSERT_EQ(basis.getID(), basis3.getID());

  ASSERT_TRUE(basis == basis2);

  ASSERT_TRUE(basis == basis3);
}

} // namespace Integrals
} // namespace Utils
} // namespace Scine