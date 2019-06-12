/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Bonds/BondOrderCollection.h>
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
  arbitraryResults.setDescription(arbitraryDescription);
}

TEST_F(AResultsTest, CanSetEnergy) {
  arbitraryResults.setEnergy(arbitraryEnergy);
}

TEST_F(AResultsTest, CanSetGradients) {
  arbitraryResults.setGradients(arbitraryGradients);
}

TEST_F(AResultsTest, CanSetHessian) {
  arbitraryResults.setHessian(arbitraryHessian);
}

TEST_F(AResultsTest, CanSetDipole) {
  arbitraryResults.setDipole(arbitraryDipole);
}

TEST_F(AResultsTest, CanSetDipoleGradient) {
  arbitraryResults.setDipoleGradient(arbitraryDipoleGradient);
}

TEST_F(AResultsTest, CanSetAODipoleMatrix) {
  arbitraryResults.setAODipoleMatrix(arbitraryDipoleAO);
}

TEST_F(AResultsTest, CanSetMODipoleMatrix) {
  arbitraryResults.setMODipoleMatrix(arbitraryDipoleMO);
}

TEST_F(AResultsTest, CanSetOneElectronMatrix) {
  arbitraryResults.setOneElectronMatrix(arbitraryOneElectronMatrix);
}

TEST_F(AResultsTest, CanSetTwoElectronMatrix) {
  arbitraryResults.setTwoElectronMatrix(arbitraryTwoElectronMatrix);
}

TEST_F(AResultsTest, CanSetBondOrders) {
  arbitraryResults.setBondOrders(arbitraryBondOrders);
}

TEST_F(AResultsTest, CanGetDescription) {
  arbitraryResults.setDescription(arbitraryDescription);
  ASSERT_EQ(arbitraryResults.getDescription(), arbitraryDescription);
}

TEST_F(AResultsTest, CanGetEnergy) {
  arbitraryResults.setEnergy(arbitraryEnergy);
  ASSERT_EQ(arbitraryResults.getEnergy(), arbitraryEnergy);
}

TEST_F(AResultsTest, CanGetGradients) {
  arbitraryResults.setGradients(arbitraryGradients);
  ASSERT_EQ(arbitraryResults.getGradients(), arbitraryGradients);
}

TEST_F(AResultsTest, CanGetHessian) {
  arbitraryResults.setHessian(arbitraryHessian);
  ASSERT_EQ(arbitraryResults.getHessian(), arbitraryHessian);
}

TEST_F(AResultsTest, CanGetDipole) {
  arbitraryResults.setDipole(arbitraryDipole);
  ASSERT_EQ(arbitraryResults.getDipole(), arbitraryDipole);
}

TEST_F(AResultsTest, CanGetDipoleGradient) {
  arbitraryResults.setDipoleGradient(arbitraryDipoleGradient);
  ASSERT_EQ(arbitraryResults.getDipoleGradient(), arbitraryDipoleGradient);
}

TEST_F(AResultsTest, CanGetAODipoleMatrix) {
  arbitraryResults.setAODipoleMatrix(arbitraryDipoleAO);
  for (int dimension = 0; dimension < 3; ++dimension) {
    ASSERT_EQ(arbitraryResults.getAODipoleMatrix()[dimension], arbitraryDipoleAO[dimension]);
  }
}

TEST_F(AResultsTest, CanGetMODipoleMatrix) {
  arbitraryResults.setMODipoleMatrix(arbitraryDipoleMO);
  for (int dimension = 0; dimension < 3; ++dimension) {
    ASSERT_EQ(arbitraryResults.getMODipoleMatrix()[dimension], arbitraryDipoleMO[dimension]);
  }
}

TEST_F(AResultsTest, CanGetOneElectronMatrix) {
  arbitraryResults.setOneElectronMatrix(arbitraryOneElectronMatrix);
  ASSERT_EQ(arbitraryResults.getOneElectronMatrix(), arbitraryOneElectronMatrix);
}

TEST_F(AResultsTest, CanGetTwoElectronMatrix) {
  arbitraryResults.setTwoElectronMatrix(arbitraryTwoElectronMatrix);
  ASSERT_EQ(arbitraryResults.getTwoElectronMatrix().restrictedMatrix(), arbitraryTwoElectronMatrix.restrictedMatrix());
}

TEST_F(AResultsTest, CanGetBondOrders) {
  arbitraryResults.setBondOrders(arbitraryBondOrders);
  ASSERT_EQ(arbitraryResults.getBondOrders(), arbitraryBondOrders);
}

TEST_F(AResultsTest, CanTakeGradients) {
  arbitraryResults.setGradients(arbitraryGradients);
  ASSERT_EQ(arbitraryResults.takeGradients(), arbitraryGradients);
  ASSERT_THROW(arbitraryResults.takeGradients(), std::exception);
}

TEST_F(AResultsTest, CanTakeHessian) {
  arbitraryResults.setHessian(arbitraryHessian);
  ASSERT_EQ(arbitraryResults.takeHessian(), arbitraryHessian);
  ASSERT_THROW(arbitraryResults.takeHessian(), std::exception);
}

TEST_F(AResultsTest, CanTakeDipoleGradients) {
  arbitraryResults.setDipoleGradient(arbitraryDipoleGradient);
  ASSERT_EQ(arbitraryResults.takeDipoleGradient(), arbitraryDipoleGradient);
  ASSERT_THROW(arbitraryResults.takeDipoleGradient(), std::exception);
}

TEST_F(AResultsTest, CanTakeAODipoleMatrix) {
  arbitraryResults.setAODipoleMatrix(arbitraryDipoleAO);
  auto aoDip = arbitraryResults.takeAODipoleMatrix();
  for (int dimension = 0; dimension < 3; ++dimension) {
    ASSERT_EQ(aoDip[dimension], arbitraryDipoleAO[dimension]);
  }
  ASSERT_THROW(arbitraryResults.takeAODipoleMatrix(), std::exception);
}

TEST_F(AResultsTest, CanTakeMODipoleMatrix) {
  arbitraryResults.setMODipoleMatrix(arbitraryDipoleMO);
  auto moDip = arbitraryResults.takeMODipoleMatrix();
  for (int dimension = 0; dimension < 3; ++dimension) {
    ASSERT_EQ(moDip[dimension], arbitraryDipoleMO[dimension]);
  }
  ASSERT_THROW(arbitraryResults.takeMODipoleMatrix(), std::exception);
}

TEST_F(AResultsTest, CanTakeOneElectronMatrix) {
  arbitraryResults.setOneElectronMatrix(arbitraryOneElectronMatrix);
  ASSERT_EQ(arbitraryResults.takeOneElectronMatrix(), arbitraryOneElectronMatrix);
  ASSERT_THROW(arbitraryResults.takeOneElectronMatrix(), std::exception);
}

TEST_F(AResultsTest, CanTakeTwoElectronMatrix) {
  arbitraryResults.setTwoElectronMatrix(arbitraryTwoElectronMatrix);
  ASSERT_EQ(arbitraryResults.takeTwoElectronMatrix().restrictedMatrix(), arbitraryTwoElectronMatrix.restrictedMatrix());
  ASSERT_THROW(arbitraryResults.takeTwoElectronMatrix(), std::exception);
}

TEST_F(AResultsTest, CanTakeBondOrders) {
  arbitraryResults.setBondOrders(arbitraryBondOrders);
  ASSERT_EQ(arbitraryResults.takeBondOrders(), arbitraryBondOrders);
  ASSERT_THROW(arbitraryResults.takeBondOrders(), std::exception);
}

TEST_F(AResultsTest, CanCopyConstructResults) {
  arbitraryResults.setEnergy(arbitraryEnergy);
  arbitraryResults.setTwoElectronMatrix(arbitraryTwoElectronMatrix);

  auto copiedResults = arbitraryResults;

  ASSERT_EQ(copiedResults.getEnergy(), arbitraryResults.getEnergy());
  ASSERT_EQ(copiedResults.getTwoElectronMatrix().restrictedMatrix(),
            arbitraryResults.getTwoElectronMatrix().restrictedMatrix());
}

TEST_F(AResultsTest, CanCopyAssignResults) {
  arbitraryResults.setEnergy(arbitraryEnergy);
  arbitraryResults.setTwoElectronMatrix(arbitraryTwoElectronMatrix);

  Results copiedResults;
  copiedResults = arbitraryResults;

  ASSERT_EQ(copiedResults.getEnergy(), arbitraryResults.getEnergy());
  ASSERT_EQ(copiedResults.getTwoElectronMatrix().restrictedMatrix(),
            arbitraryResults.getTwoElectronMatrix().restrictedMatrix());
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
