/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/CalculatorBasics/PropertyList.h>
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Utils {
namespace Tests {

class APropertyListTest : public Test {
 public:
  PropertyList aPropertyList;
};

TEST_F(APropertyListTest, CanAddEnergy) {
  aPropertyList.addProperty(Property::Energy);
  ASSERT_TRUE(aPropertyList.containsSubSet(Property::Energy));
}

TEST_F(APropertyListTest, CanAddGradient) {
  aPropertyList.addProperty(Property::Gradients);
  ASSERT_TRUE(aPropertyList.containsSubSet(Property::Gradients));
}

TEST_F(APropertyListTest, CanAddHessian) {
  aPropertyList.addProperty(Property::Hessian);
}

TEST_F(APropertyListTest, CanAddPartialHessian) {
  aPropertyList.addProperty(Property::PartialHessian);
}

TEST_F(APropertyListTest, CanAddMultipleProperties) {
  auto a = PropertyList(Property::Energy);
  a.addProperty(Property::Gradients);
  a.addProperty(Property::Hessian);
  auto b = PropertyList();
  b.addProperties(a);
  ASSERT_TRUE(a.containsSubSet(Property::Energy) && a.containsSubSet(Property::Gradients) &&
              a.containsSubSet(Property::Hessian));
  ASSERT_TRUE(b.containsSubSet(Property::Energy) && b.containsSubSet(Property::Gradients) &&
              b.containsSubSet(Property::Hessian));
}

TEST_F(APropertyListTest, CanGetIntersection) {
  auto a = PropertyList(Property::Energy);
  a.addProperty(Property::Gradients);
  a.addProperty(Property::Hessian);
  auto b = PropertyList(Property::Energy | Property::Hessian | Property::AtomicCharges);
  auto c = a.intersection(b);
  ASSERT_TRUE(a.containsSubSet(Property::Energy) && a.containsSubSet(Property::Gradients) &&
              a.containsSubSet(Property::Hessian));
  ASSERT_TRUE(b.containsSubSet(Property::Energy) && !b.containsSubSet(Property::Gradients) &&
              b.containsSubSet(Property::Hessian) && b.containsSubSet(Property::AtomicCharges));
  ASSERT_TRUE(c.containsSubSet(Property::Energy) && !c.containsSubSet(Property::Gradients) &&
              c.containsSubSet(Property::Hessian) && !c.containsSubSet(Property::AtomicCharges));
}

TEST_F(APropertyListTest, CanAddMolecularDipole) {
  aPropertyList.addProperty(Property::Dipole);
}

TEST_F(APropertyListTest, CanAddDipoleMatrices) {
  aPropertyList.addProperty(Property::DipoleMatrixAO | Property::DipoleMatrixMO);
}

TEST_F(APropertyListTest, CanAddElectronMatrices) {
  aPropertyList.addProperty(Property::OneElectronMatrix | Property::TwoElectronMatrix);
}

TEST_F(APropertyListTest, CanAddEverything) {
  aPropertyList.addProperty(
      Property::Energy | Property::Gradients | Property::Hessian | Property::PartialHessian | Property::AtomicHessians |
      Property::Dipole | Property::DipoleGradient | Property::DipoleMatrixAO | Property::DipoleMatrixMO |
      Property::DensityMatrix | Property::OneElectronMatrix | Property::TwoElectronMatrix | Property::OverlapMatrix |
      Property::CoefficientMatrix | Property::OrbitalEnergies | Property::ElectronicOccupation | Property::Thermochemistry |
      Property::ExcitedStates | Property::AOtoAtomMapping | Property::AtomicCharges | Property::BondOrderMatrix |
      Property::Description | Property::SuccessfulCalculation | Property::ProgramName | Property::PointChargesGradients |
      Property::MoessbauerParameter | Property::AtomicGtos | Property::GridOccupation | Property::StressTensor);
}

TEST_F(APropertyListTest, CanRemoveAProperty) {
  PropertyList newPropertyList;
  newPropertyList.addProperty(Property::Energy | Property::Hessian);
  newPropertyList.removeProperty(Property::Hessian);
  ASSERT_FALSE(newPropertyList.containsSubSet(Utils::Property::Hessian));
  ASSERT_TRUE(newPropertyList.containsSubSet(Utils::Property::Energy));
  newPropertyList.addProperty(Property::Hessian);
  newPropertyList.addProperty(Property::Gradients);
  ASSERT_TRUE(newPropertyList.containsSubSet(Utils::Property::Hessian));
  ASSERT_TRUE(newPropertyList.containsSubSet(Utils::Property::Gradients));
  newPropertyList.removeProperties(Property::Hessian | Property::Energy);
  ASSERT_FALSE(newPropertyList.containsSubSet(Utils::Property::Hessian));
  ASSERT_FALSE(newPropertyList.containsSubSet(Utils::Property::Energy));
  ASSERT_TRUE(newPropertyList.containsSubSet(Utils::Property::Gradients));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
