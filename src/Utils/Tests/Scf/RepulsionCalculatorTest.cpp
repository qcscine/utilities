/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Utils/Math/FullSecondDerivativeCollection.h>
#include <Utils/Scf/MethodInterfaces/RepulsionCalculator.h>
#include <Utils/Typenames.h>
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Utils {

class ACoreCoreRepulsionCalculatorTest : public Test {
 public:
  ElementTypeCollection arbitraryElements;
  PositionCollection arbitraryPositions;
  std::unique_ptr<RepulsionCalculator> repulsionCalculator;
  HessianMatrix calculateFromEnergyDifferences(double delta);

 private:
  double hessianElementSameFromEnergy(int i, const PositionCollection& referencePositions, double delta);
  double hessianElementDifferentFromEnergy(int i, int j, const PositionCollection& referencePositions, double delta);

  void SetUp() final {
    arbitraryElements = {ElementType::C, ElementType::H, ElementType::H, ElementType::H, ElementType::H};
    arbitraryPositions = PositionCollection::Random(arbitraryElements.size(), 3);

    repulsionCalculator = std::make_unique<RepulsionCalculator>(arbitraryElements, arbitraryPositions);
  }
};

HessianMatrix ACoreCoreRepulsionCalculatorTest::calculateFromEnergyDifferences(double delta) {
  auto positions = arbitraryPositions;

  auto nDimensions = 3 * positions.rows();
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(nDimensions, nDimensions);

  for (int i = 0; i < nDimensions; ++i) {
    m(i, i) = hessianElementSameFromEnergy(i, positions, delta);
    for (int j = 0; j < i; ++j) {
      m(i, j) = hessianElementDifferentFromEnergy(i, j, positions, delta);
      m(j, i) = m(i, j);
    }
  }

  return HessianMatrix{m};
}

double ACoreCoreRepulsionCalculatorTest::hessianElementSameFromEnergy(int i, const PositionCollection& referencePositions,
                                                                      double delta) {
  auto atomI = i / 3;
  auto compI = i % 3;

  arbitraryPositions = referencePositions;
  repulsionCalculator->calculateRepulsion(DerivativeOrder::Zero);
  auto E = repulsionCalculator->getRepulsionEnergy();

  arbitraryPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) - delta;
  repulsionCalculator->calculateRepulsion(DerivativeOrder::Zero);
  auto Em = repulsionCalculator->getRepulsionEnergy();

  arbitraryPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) + delta;
  repulsionCalculator->calculateRepulsion(DerivativeOrder::Zero);
  auto Ep = repulsionCalculator->getRepulsionEnergy();

  return (Ep - 2 * E + Em) / (delta * delta);
}

double ACoreCoreRepulsionCalculatorTest::hessianElementDifferentFromEnergy(int i, int j,
                                                                           const PositionCollection& referencePositions,
                                                                           double delta) {
  auto D = delta / 2;

  auto atomI = i / 3;
  auto compI = i % 3;
  auto atomJ = j / 3;
  auto compJ = j % 3;

  arbitraryPositions = referencePositions;

  arbitraryPositions.row(atomI)[compI] = referencePositions.row(atomI)(compI) + D;
  arbitraryPositions.row(atomJ)(compJ) = referencePositions.row(atomJ)(compJ) + D;
  repulsionCalculator->calculateRepulsion(DerivativeOrder::Zero);
  auto Epp = repulsionCalculator->getRepulsionEnergy();

  arbitraryPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) - D;
  arbitraryPositions.row(atomJ)(compJ) = referencePositions.row(atomJ)(compJ) + D;
  repulsionCalculator->calculateRepulsion(DerivativeOrder::Zero);
  auto Emp = repulsionCalculator->getRepulsionEnergy();

  arbitraryPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) + D;
  arbitraryPositions.row(atomJ)(compJ) = referencePositions.row(atomJ)(compJ) - D;
  repulsionCalculator->calculateRepulsion(DerivativeOrder::Zero);
  auto Epm = repulsionCalculator->getRepulsionEnergy();

  arbitraryPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) - D;
  arbitraryPositions.row(atomJ)(compJ) = referencePositions.row(atomJ)(compJ) - D;
  repulsionCalculator->calculateRepulsion(DerivativeOrder::Zero);
  auto Emm = repulsionCalculator->getRepulsionEnergy();

  return (Epp - Epm - Emp + Emm) / (delta * delta);
}
TEST_F(ACoreCoreRepulsionCalculatorTest, CanCalculateRepulsionEnergy) {
  double repulsionEnergy = 0;
  for (int i = 0; i < arbitraryPositions.rows(); ++i) {
    for (int j = i + 1; j < arbitraryPositions.rows(); ++j) {
      double repulsionConstant = ElementInfo::Z(arbitraryElements[i]) * ElementInfo::Z(arbitraryElements[j]);
      double distance = (arbitraryPositions.row(i) - arbitraryPositions.row(j)).norm();
      repulsionEnergy += repulsionConstant / distance;
    }
  }

  repulsionCalculator->calculateRepulsion(DerivativeOrder::Zero);
  ASSERT_DOUBLE_EQ(repulsionCalculator->getRepulsionEnergy(), repulsionEnergy);
}

TEST_F(ACoreCoreRepulsionCalculatorTest, CanCalculateRepulsionGradient) {
  repulsionCalculator->calculateRepulsion(DerivativeOrder::One);
  GradientCollection gradients = GradientCollection::Zero(arbitraryElements.size(), 3);
  repulsionCalculator->addRepulsionDerivatives(gradients);
  for (int i = 0; i < arbitraryPositions.rows(); ++i) {
    for (int j = 0; j < arbitraryPositions.cols(); ++j) {
      const double position = arbitraryPositions(i, j);
      arbitraryPositions(i, j) += 1e-5;
      repulsionCalculator->calculateRepulsion(DerivativeOrder::Zero);
      const double energyPlus = repulsionCalculator->getRepulsionEnergy();
      arbitraryPositions(i, j) = position - 1e-5;
      repulsionCalculator->calculateRepulsion(DerivativeOrder::Zero);
      const double energyMinus = repulsionCalculator->getRepulsionEnergy();
      arbitraryPositions(i, j) = position;
      EXPECT_THAT((energyPlus - energyMinus) / 2e-5, DoubleNear(gradients(i, j), 1e-5));
    }
  }
}

TEST_F(ACoreCoreRepulsionCalculatorTest, CanCalculateRepulsionHessian) {
  repulsionCalculator->calculateRepulsion(DerivativeOrder::Two);
  FullSecondDerivativeCollection secondDerivatives(arbitraryElements.size());
  secondDerivatives.setZero();
  repulsionCalculator->addRepulsionDerivatives(secondDerivatives);
  HessianMatrix hessian = secondDerivatives.getHessianMatrix();

  using MapType = Eigen::Map<const Eigen::VectorXd>;

  HessianMatrix numericalHessian = HessianMatrix::Zero(arbitraryPositions.size(), arbitraryPositions.size());
  // HessianMatrix numericalHessian = calculateFromEnergyDifferences(1e-6);

  HessianMatrix::Zero(3 * arbitraryElements.size(), 3 * arbitraryElements.size());
  for (unsigned hessianIndex = 0; hessianIndex < 3 * arbitraryElements.size(); ++hessianIndex) {
    const int i = hessianIndex / 3;
    const int j = hessianIndex % 3;
    const double position = arbitraryPositions(i, j);
    arbitraryPositions(i, j) += 1e-6;
    repulsionCalculator->calculateRepulsion(DerivativeOrder::One);
    GradientCollection gradPlus = GradientCollection::Zero(arbitraryElements.size(), 3);
    repulsionCalculator->addRepulsionDerivatives(gradPlus);
    const Eigen::VectorXd mappedGradientPlus = MapType(gradPlus.data(), gradPlus.size());
    arbitraryPositions(i, j) = position - 1e-6;
    repulsionCalculator->calculateRepulsion(DerivativeOrder::One);
    GradientCollection gradMinus = GradientCollection::Zero(arbitraryElements.size(), 3);
    repulsionCalculator->addRepulsionDerivatives(gradMinus);
    const Eigen::VectorXd mappedGradientMinus = MapType(gradMinus.data(), gradMinus.size());
    arbitraryPositions(i, j) = position;
    numericalHessian.row(hessianIndex) = (mappedGradientPlus - mappedGradientMinus) / 2e-6;
  }

  for (int i = 0; i < hessian.size(); ++i) {
    EXPECT_THAT(hessian(i), DoubleNear(numericalHessian(i), 1e-7));
  }
}
} // namespace Utils
} // namespace Scine
