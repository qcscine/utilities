/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Constants.h>
#include <Utils/Geometry/ElementInfo.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

static_assert(ElementInfo::Z(ElementType::H) == 1, "H A isn't right");
static_assert(ElementInfo::A(ElementType::H) == 0, "H A isn't right");
static_assert(ElementInfo::Z(ElementType::H1) == 1, "H1 A isn't right");
static_assert(ElementInfo::A(ElementType::H1) == 1, "H1 A isn't right");
static_assert(ElementInfo::Z(ElementType::Be) == 4, "Be Z isn't right");
static_assert(ElementInfo::A(ElementType::Be) == 9, "Be A isn't right");
static_assert(ElementInfo::Z(ElementType::V) == 23, "V Z isn't right");
static_assert(ElementInfo::A(ElementType::V) == 0, "V A isn't right");
static_assert(ElementInfo::Z(ElementType::V50) == 23, "V50 Z isn't right");
static_assert(ElementInfo::A(ElementType::V50) == 50, "V50 A isn't right");

TEST(ElementInfoTest, CanReadMassesInParallel) {
  double massSum = 0.;
  double mass100Hydrogen = 1.0079 * 100;

#pragma omp parallel for reduction(+ : massSum)
  for (int i = 0; i < 100; ++i) {
    massSum += ElementInfo::mass(ElementType::H);
  }
  ASSERT_THAT(massSum, DoubleNear(mass100Hydrogen, 1e-11));
}

TEST(ElementInfoTest, CanConvertSymbolStringToElementType) {
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("h"), Eq(ElementType::H));
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("H"), Eq(ElementType::H));
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("1H"), Eq(ElementType::H1));
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("H1"), Eq(ElementType::H1));
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("D"), Eq(ElementType::D));
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("2H"), Eq(ElementType::D));
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("T"), Eq(ElementType::T));
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("3H"), Eq(ElementType::T));
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("C"), Eq(ElementType::C));
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("C12"), Eq(ElementType::C12));
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("12C"), Eq(ElementType::C12));
  ASSERT_THAT(ElementInfo::elementTypeForSymbol("Rn"), Eq(ElementType::Rn));
}

TEST(ElementInfoTest, CanConvertElementTypeToSymbol) {
  ASSERT_THAT(ElementInfo::symbol(ElementType::H), Eq(std::string("H")));
  ASSERT_THAT(ElementInfo::symbol(ElementType::H1), Eq(std::string("H")));
  ASSERT_THAT(ElementInfo::symbol(ElementType::C), Eq(std::string("C")));
  ASSERT_THAT(ElementInfo::symbol(ElementType::C12), Eq(std::string("C")));
  ASSERT_THAT(ElementInfo::symbol(ElementType::Rn), Eq(std::string("Rn")));
}

TEST(ElementInfoTest, CanGetVdwRadiusForElementInAtomicUnits) {
  double vdwInPicometers = 210.0;
  double vdwInAngstrom = vdwInPicometers / 100.0;
  double vdwInAtomicUnits = vdwInAngstrom * Constants::bohr_per_angstrom;
  ASSERT_THAT(ElementInfo::vdwRadius(ElementType::Si), Eq(vdwInAtomicUnits));
}

TEST(ElementInfoTest, CanGetCovalentRadiusForElementInAtomicUnits) {
  double covalentInPicometers = 134.0;
  double covalentInAngstrom = covalentInPicometers / 100.0;
  double covalentInAtomicUnits = covalentInAngstrom * Constants::bohr_per_angstrom;
  ASSERT_THAT(ElementInfo::covalentRadius(ElementType::Hs), Eq(covalentInAtomicUnits));
}

TEST(ElementInfoTest, IsotopeMasses) {
  ASSERT_THAT(ElementInfo::mass(ElementType::V), Eq(50.942));
  ASSERT_THAT(ElementInfo::mass(ElementType::V50), Eq(49.94715601));
}

TEST(ElementInfoTest, IsotopeBase) {
  // Identity and base tests for non-monoisotopic elements
  ASSERT_THAT(ElementInfo::base(ElementType::H), Eq(ElementType::H));
  ASSERT_THAT(ElementInfo::base(ElementType::H1), Eq(ElementType::H));
  ASSERT_THAT(ElementInfo::base(ElementType::D), Eq(ElementType::H));
  ASSERT_THAT(ElementInfo::base(ElementType::V50), Eq(ElementType::V));
  ASSERT_THAT(ElementInfo::base(ElementType::V), Eq(ElementType::V));
  // Identity tests for monoisotopic element
  ASSERT_THAT(ElementInfo::base(ElementType::Mn55), Eq(ElementType::Mn));
  ASSERT_THAT(ElementInfo::base(ElementType::Mn), Eq(ElementType::Mn));
}

TEST(ElementInfoTest, ElementIsotopes) {
  const std::vector<ElementType> hydrogenIsotopes{ElementType::H1, ElementType::D, ElementType::T};

  // Isotopes are unordered, so for comparison we need to sort them
  auto sortedIsotopes = [](const ElementType e) {
    auto isotopes = ElementInfo::isotopes(e);
    std::sort(std::begin(isotopes), std::end(isotopes));
    return isotopes;
  };

  ASSERT_THAT(sortedIsotopes(ElementType::H), Eq(hydrogenIsotopes));
  ASSERT_THAT(sortedIsotopes(ElementType::H1), Eq(hydrogenIsotopes));
  // Monoisotopic element
  ASSERT_THAT(ElementInfo::isotopes(ElementType::Co59), Eq(std::vector<ElementType>(1, ElementType::Co)));
}

TEST(ElementInfoTest, ComposeElements) {
  ASSERT_THAT(ElementInfo::element(3), Eq(ElementType::Li));
  ASSERT_THAT(ElementInfo::element(4), Eq(ElementType::Be));
  ASSERT_THAT(ElementInfo::element(112), Eq(ElementType::Cn));

  ASSERT_THAT(ElementInfo::isotope(98, 251), Eq(ElementType::Cf251));
  ASSERT_THAT(ElementInfo::isotope(1, 1), Eq(ElementType::H1));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
