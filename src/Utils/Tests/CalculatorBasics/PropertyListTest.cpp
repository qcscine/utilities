/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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
}

TEST_F(APropertyListTest, CanAddGradient) {
  aPropertyList.addProperty(Property::Gradients);
}

TEST_F(APropertyListTest, CanAddHessian) {
  aPropertyList.addProperty(Property::Hessian);
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
  aPropertyList.addProperty(Property::Energy | Property::Gradients | Property::Hessian | Property::Dipole |
                            Property::DipoleMatrixAO | Property::DipoleMatrixMO | Property::OneElectronMatrix |
                            Property::TwoElectronMatrix);
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
