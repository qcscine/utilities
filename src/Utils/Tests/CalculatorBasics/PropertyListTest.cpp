/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
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
  aPropertyList.addProperty(
      Property::Energy | Property::Gradients | Property::Hessian | Property::AtomicHessians | Property::Dipole |
      Property::DipoleGradient | Property::DipoleMatrixAO | Property::DipoleMatrixMO | Property::DensityMatrix |
      Property::OneElectronMatrix | Property::TwoElectronMatrix | Property::OverlapMatrix |
      Property::CoefficientMatrix | Property::OrbitalEnergies | Property::ElectronicOccupation |
      Property::Thermochemistry | Property::ExcitedStates | Property::AOtoAtomMapping | Property::AtomicCharges |
      Property::BondOrderMatrix | Property::Description | Property::SuccessfulCalculation | Property::ProgramName |
      Property::PointChargesGradients | Property::AtomicGtos | Property::GridOccupation | Property::StressTensor);
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
