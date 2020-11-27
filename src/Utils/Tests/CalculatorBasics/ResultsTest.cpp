/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Typenames.h>
#include <gmock/gmock.h>
#include <Eigen/Core>

using namespace testing;

namespace Scine {
namespace Utils {
namespace Tests {

class AResultsTest : public Test {
 public:
  Results arbitraryResults;
  std::string arbitraryDescription = "";
  double arbitraryEnergy = -145.34;
  GradientCollection arbitraryGradients;
  HessianMatrix arbitraryHessian;
  Dipole arbitraryDipole;
  DipoleGradient arbitraryDipoleGradient;
  DipoleMatrix arbitraryDipoleAO, arbitraryDipoleMO;
  Eigen::MatrixXd arbitraryOneElectronMatrix;
  SpinAdaptedMatrix arbitraryTwoElectronMatrix;
  BondOrderCollection arbitraryBondOrders;

  int arbitraryDimension = 10;
  int arbitraryAtomNumber = 4;

 private:
  void SetUp() final {
    arbitraryGradients = GradientCollection::Random(arbitraryAtomNumber, 3);
    arbitraryHessian = HessianMatrix::Random(arbitraryAtomNumber * 3, arbitraryAtomNumber * 3);
    arbitraryDipole = Dipole::Random(3);
    arbitraryDipoleGradient = DipoleGradient::Random(arbitraryAtomNumber * 3, 3);
    for (int dimension = 0; dimension < 3; ++dimension) {
      arbitraryDipoleAO[dimension] = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
      arbitraryDipoleMO[dimension] = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    }
    arbitraryOneElectronMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    arbitraryTwoElectronMatrix =
        SpinAdaptedMatrix::createRestricted(Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension));

    arbitraryBondOrders.resize(arbitraryDimension);
    {
      // Initialize sparse bond order matrix using a dense random matrix
      Eigen::MatrixXd denseRandom = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);

      for (Eigen::Index i = 0; i < arbitraryDimension; ++i) {
        for (Eigen::Index j = i + 1; j < arbitraryDimension; ++j) {
          double randomValue = denseRandom(i, j);
          // Use 25% of the values for sparseness (value range is [-1, 1])
          if (randomValue > 0.5) {
            arbitraryBondOrders.setOrder(i, j, randomValue);
          }
        }
      }
    }
  }
};

TEST_F(AResultsTest, CanCreateDefaultResult) {
  Results result;
}

TEST_F(AResultsTest, CanSetDescription) {
  arbitraryResults.set<Property::Description>(arbitraryDescription);
}

TEST_F(AResultsTest, CanSetEnergy) {
  arbitraryResults.set<Property::Energy>(arbitraryEnergy);
}

TEST_F(AResultsTest, CanSetGradients) {
  arbitraryResults.set<Property::Gradients>(arbitraryGradients);
}

TEST_F(AResultsTest, CanSetHessian) {
  arbitraryResults.set<Property::Hessian>(arbitraryHessian);
}

TEST_F(AResultsTest, CanSetDipole) {
  arbitraryResults.set<Property::Dipole>(arbitraryDipole);
}

TEST_F(AResultsTest, CanSetDipoleGradient) {
  arbitraryResults.set<Property::DipoleGradient>(arbitraryDipoleGradient);
}

TEST_F(AResultsTest, CanSetAODipoleMatrix) {
  arbitraryResults.set<Property::DipoleMatrixAO>(arbitraryDipoleAO);
}

TEST_F(AResultsTest, CanSetMODipoleMatrix) {
  arbitraryResults.set<Property::DipoleMatrixMO>(arbitraryDipoleMO);
}

TEST_F(AResultsTest, CanSetOneElectronMatrix) {
  arbitraryResults.set<Property::OneElectronMatrix>(arbitraryOneElectronMatrix);
}

TEST_F(AResultsTest, CanSetTwoElectronMatrix) {
  arbitraryResults.set<Property::TwoElectronMatrix>(arbitraryTwoElectronMatrix);
}

TEST_F(AResultsTest, CanSetBondOrders) {
  arbitraryResults.set<Property::BondOrderMatrix>(arbitraryBondOrders);
}

TEST_F(AResultsTest, CanGetDescription) {
  arbitraryResults.set<Property::Description>(arbitraryDescription);
  ASSERT_EQ(arbitraryResults.get<Property::Description>(), arbitraryDescription);
}

TEST_F(AResultsTest, CanGetEnergy) {
  arbitraryResults.set<Property::Energy>(arbitraryEnergy);
  ASSERT_EQ(arbitraryResults.get<Property::Energy>(), arbitraryEnergy);
}

TEST_F(AResultsTest, CanGetGradients) {
  arbitraryResults.set<Property::Gradients>(arbitraryGradients);
  ASSERT_EQ(arbitraryResults.get<Property::Gradients>(), arbitraryGradients);
}

TEST_F(AResultsTest, CanGetHessian) {
  arbitraryResults.set<Property::Hessian>(arbitraryHessian);
  ASSERT_EQ(arbitraryResults.get<Property::Hessian>(), arbitraryHessian);
}

TEST_F(AResultsTest, CanGetDipole) {
  arbitraryResults.set<Property::Dipole>(arbitraryDipole);
  ASSERT_EQ(arbitraryResults.get<Property::Dipole>(), arbitraryDipole);
}

TEST_F(AResultsTest, CanGetDipoleGradient) {
  arbitraryResults.set<Property::DipoleGradient>(arbitraryDipoleGradient);
  ASSERT_EQ(arbitraryResults.get<Property::DipoleGradient>(), arbitraryDipoleGradient);
}

TEST_F(AResultsTest, CanGetAODipoleMatrix) {
  arbitraryResults.set<Property::DipoleMatrixAO>(arbitraryDipoleAO);
  for (int dimension = 0; dimension < 3; ++dimension) {
    ASSERT_EQ(arbitraryResults.get<Property::DipoleMatrixAO>()[dimension], arbitraryDipoleAO[dimension]);
  }
}

TEST_F(AResultsTest, CanGetMODipoleMatrix) {
  arbitraryResults.set<Property::DipoleMatrixMO>(arbitraryDipoleMO);
  for (int dimension = 0; dimension < 3; ++dimension) {
    ASSERT_EQ(arbitraryResults.get<Property::DipoleMatrixMO>()[dimension], arbitraryDipoleMO[dimension]);
  }
}

TEST_F(AResultsTest, CanGetOneElectronMatrix) {
  arbitraryResults.set<Property::OneElectronMatrix>(arbitraryOneElectronMatrix);
  ASSERT_EQ(arbitraryResults.get<Property::OneElectronMatrix>(), arbitraryOneElectronMatrix);
}

TEST_F(AResultsTest, CanGetTwoElectronMatrix) {
  arbitraryResults.set<Property::TwoElectronMatrix>(arbitraryTwoElectronMatrix);
  ASSERT_EQ(arbitraryResults.get<Property::TwoElectronMatrix>().restrictedMatrix(),
            arbitraryTwoElectronMatrix.restrictedMatrix());
}

TEST_F(AResultsTest, CanGetBondOrders) {
  arbitraryResults.set<Property::BondOrderMatrix>(arbitraryBondOrders);
  ASSERT_EQ(arbitraryResults.get<Property::BondOrderMatrix>(), arbitraryBondOrders);
}

TEST_F(AResultsTest, CanTakeGradients) {
  arbitraryResults.set<Property::Gradients>(arbitraryGradients);
  ASSERT_EQ(arbitraryResults.take<Property::Gradients>(), arbitraryGradients);
  ASSERT_THROW(arbitraryResults.take<Property::Gradients>(), std::exception);
}

TEST_F(AResultsTest, CanTakeHessian) {
  arbitraryResults.set<Property::Hessian>(arbitraryHessian);
  ASSERT_EQ(arbitraryResults.take<Property::Hessian>(), arbitraryHessian);
  ASSERT_THROW(arbitraryResults.take<Property::Hessian>(), std::exception);
}

TEST_F(AResultsTest, CanTakeDipoleGradients) {
  arbitraryResults.set<Property::DipoleGradient>(arbitraryDipoleGradient);
  ASSERT_EQ(arbitraryResults.take<Property::DipoleGradient>(), arbitraryDipoleGradient);
  ASSERT_THROW(arbitraryResults.take<Property::DipoleGradient>(), std::exception);
}

TEST_F(AResultsTest, CanTakeAODipoleMatrix) {
  arbitraryResults.set<Property::DipoleMatrixAO>(arbitraryDipoleAO);
  auto aoDip = arbitraryResults.take<Property::DipoleMatrixAO>();
  for (int dimension = 0; dimension < 3; ++dimension) {
    ASSERT_EQ(aoDip[dimension], arbitraryDipoleAO[dimension]);
  }
  ASSERT_THROW(arbitraryResults.take<Property::DipoleMatrixAO>(), std::exception);
}

TEST_F(AResultsTest, CanTakeMODipoleMatrix) {
  arbitraryResults.set<Property::DipoleMatrixMO>(arbitraryDipoleMO);
  auto moDip = arbitraryResults.take<Property::DipoleMatrixMO>();
  for (int dimension = 0; dimension < 3; ++dimension) {
    ASSERT_EQ(moDip[dimension], arbitraryDipoleMO[dimension]);
  }
  ASSERT_THROW(arbitraryResults.take<Property::DipoleMatrixMO>(), std::exception);
}

TEST_F(AResultsTest, CanTakeOneElectronMatrix) {
  arbitraryResults.set<Property::OneElectronMatrix>(arbitraryOneElectronMatrix);
  ASSERT_EQ(arbitraryResults.take<Property::OneElectronMatrix>(), arbitraryOneElectronMatrix);
  ASSERT_THROW(arbitraryResults.take<Property::OneElectronMatrix>(), std::exception);
}

TEST_F(AResultsTest, CanTakeTwoElectronMatrix) {
  arbitraryResults.set<Property::TwoElectronMatrix>(arbitraryTwoElectronMatrix);
  ASSERT_EQ(arbitraryResults.take<Property::TwoElectronMatrix>().restrictedMatrix(),
            arbitraryTwoElectronMatrix.restrictedMatrix());
  ASSERT_THROW(arbitraryResults.take<Property::TwoElectronMatrix>(), std::exception);
}

TEST_F(AResultsTest, CanTakeBondOrders) {
  arbitraryResults.set<Property::BondOrderMatrix>(arbitraryBondOrders);
  ASSERT_EQ(arbitraryResults.take<Property::BondOrderMatrix>(), arbitraryBondOrders);
  ASSERT_THROW(arbitraryResults.take<Property::BondOrderMatrix>(), std::exception);
}

TEST_F(AResultsTest, CanCopyConstructResults) {
  arbitraryResults.set<Property::Energy>(arbitraryEnergy);
  arbitraryResults.set<Property::TwoElectronMatrix>(arbitraryTwoElectronMatrix);

  auto copiedResults = arbitraryResults;

  ASSERT_EQ(copiedResults.get<Property::Energy>(), arbitraryResults.get<Property::Energy>());
  ASSERT_EQ(copiedResults.get<Property::TwoElectronMatrix>().restrictedMatrix(),
            arbitraryResults.get<Property::TwoElectronMatrix>().restrictedMatrix());
}

TEST_F(AResultsTest, CanCopyAssignResults) {
  arbitraryResults.set<Property::Energy>(arbitraryEnergy);
  arbitraryResults.set<Property::TwoElectronMatrix>(arbitraryTwoElectronMatrix);

  Results copiedResults;
  copiedResults = arbitraryResults;

  ASSERT_EQ(copiedResults.get<Property::Energy>(), arbitraryResults.get<Property::Energy>());
  ASSERT_EQ(copiedResults.get<Property::TwoElectronMatrix>().restrictedMatrix(),
            arbitraryResults.get<Property::TwoElectronMatrix>().restrictedMatrix());
}

TEST_F(AResultsTest, CanReturnAllContainedPropertiesList) {
  arbitraryResults.set<Property::Energy>(arbitraryEnergy);
  arbitraryResults.set<Property::TwoElectronMatrix>(arbitraryTwoElectronMatrix);

  PropertyList samePropertyList;
  samePropertyList.addProperty(Property::Energy);
  samePropertyList.addProperty(Property::TwoElectronMatrix);

  PropertyList differentPropertyList;
  differentPropertyList.addProperty(Property::Energy);
  differentPropertyList.addProperty(Property::Gradients);

  PropertyList containedPropertyList = arbitraryResults.allContainedProperties();

  ASSERT_TRUE(containedPropertyList.containsSubSet(samePropertyList) && samePropertyList.containsSubSet(containedPropertyList));
  ASSERT_FALSE(containedPropertyList.containsSubSet(differentPropertyList) &&
               differentPropertyList.containsSubSet(containedPropertyList));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
