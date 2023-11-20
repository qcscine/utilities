/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/Log.h>
#include <Utils/DataStructures/DipoleMatrix.h>
#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>
#include <Utils/Scf/MethodInterfaces/AdditiveElectronicContribution.h>
#include <Utils/Scf/MethodInterfaces/DensityMatrixGuessCalculator.h>
#include <Utils/Scf/MethodInterfaces/ElectronicContributionCalculator.h>
#include <Utils/Scf/MethodInterfaces/LcaoMethod.h>
#include <Utils/Scf/MethodInterfaces/OverlapCalculator.h>
#include <Utils/Scf/MethodInterfaces/RepulsionCalculator.h>
#include <Utils/Scf/MethodInterfaces/ScfMethod.h>
#include <Utils/Scf/MethodInterfaces/StructureDependentInitializer.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

// Line of electrons
struct TestOverlapCalculator final : public OverlapCalculator {
  ~TestOverlapCalculator() final = default;
  void calculateOverlap(Utils::DerivativeOrder /*highestRequiredOrder*/) final {
    Eigen::MatrixXd ovlp(10, 10);
    for (int i = 0; i < 10; ++i) {
      for (int j = 0; j < 10; ++j) {
        int distance = std::abs(i - j);
        ovlp(i, j) = std::exp(-1.0 * distance);
      }
    }
    returnVal_.setBaseMatrix(ovlp);
  }
  const MatrixWithDerivatives& getOverlap() const final {
    return returnVal_;
  }
  void reset() final{};

 private:
  MatrixWithDerivatives returnVal_;
};

struct TestElectronicPart final : public ElectronicContributionCalculator {
  TestElectronicPart(const DensityMatrix& density, const AtomsOrbitalsIndexes& aoIndex)
    : denmat_(density), aoIndex_(aoIndex) {
  }
  ~TestElectronicPart() final = default;
  /*! Reinitialize after a change in the elements. */
  void initialize() final{};

  /** This function will be called only once per single-point calculation. */
  void calculateDensityIndependentPart(Utils::DerivativeOrder /*order*/) final {
    TestOverlapCalculator ovlpCalc;
    ovlpCalc.calculateOverlap(Utils::DerivativeOrder::Zero);
    oneElMat = -2.5 * ovlpCalc.getOverlap().getMatrixXd();
  }
  void calculateDensityDependentPart(Utils::DerivativeOrder /*order*/) final {
    // Charges = S*P
    std::vector<double> mullikenCharges(10), coreCharges(10);
    std::fill(coreCharges.begin(), coreCharges.end(), 1.0);
    TestOverlapCalculator ovlpCalc;
    ovlpCalc.calculateOverlap(Utils::DerivativeOrder::Zero);

    LcaoUtils::calculateMullikenAtomicCharges(mullikenCharges, coreCharges, denmat_,
                                              ovlpCalc.getOverlap().getMatrixXd(), aoIndex_);

    twoElMat = Eigen::MatrixXd::Zero(10, 10);
    for (int i = 0; i < 10; ++i) {
      twoElMat(i, i) = 2.0;
      for (int j = 0; j < i; ++j) {
        twoElMat(i, j) = (mullikenCharges[i] - 1.0) * (mullikenCharges[j] - 1.0) / double(std::abs(i - j));
      }
    }
    if (fieldContribution_) {
      fieldContribution_->calculate(denmat_, DerivativeOrder::Zero);
    }
    twoElMat = twoElMat.selfadjointView<Eigen::Lower>();
  }
  void finalize(Utils::DerivativeOrder order) final {
    calculateDensityDependentPart(order);
  }
  void addDensityDependentElectronicContribution(std::shared_ptr<AdditiveElectronicContribution> contribution) final {
    fieldContribution_ = std::move(contribution);
  }
  void addDensityIndependentElectronicContribution(std::shared_ptr<AdditiveElectronicContribution> /*contribution*/) final {
  }

  SpinAdaptedMatrix getMatrix() const final {
    Eigen::MatrixXd fock = oneElMat + twoElMat;
    if (fieldContribution_) {
      fock += fieldContribution_->getElectronicContribution().restrictedMatrix();
    }
    return SpinAdaptedMatrix::createRestricted(std::move(fock));
  }

  double calculateElectronicEnergy() const final {
    auto nAOs = denmat_.restrictedMatrix().rows();
    double electronicEnergy = 0.0;
    for (unsigned int i = 0; i < nAOs; i++) {
      electronicEnergy += 0.5 * denmat_.restricted(i, i) * (twoElMat(i, i) + 2 * oneElMat(i, i));
      for (unsigned int j = 0; j < i; j++) {
        electronicEnergy += denmat_.restricted(i, j) * (twoElMat(i, j) + 2 * oneElMat(i, j));
      }
    }
    if (fieldContribution_) {
      electronicEnergy += fieldContribution_->getElectronicEnergyContribution();
    }
    return electronicEnergy;
  }
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::Derivative::First>& /*d*/) const final{};
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::Derivative::SecondAtomic>& /*d*/) const final{};
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::Derivative::SecondFull>& /*d*/) const final{};

 private:
  const DensityMatrix& denmat_;
  const AtomsOrbitalsIndexes& aoIndex_;
  Eigen::MatrixXd oneElMat;
  Eigen::MatrixXd twoElMat;
  std::shared_ptr<AdditiveElectronicContribution> fieldContribution_;
}; // namespace Tests

struct TestDenMatGuessCalc final : public DensityMatrixGuessCalculator {
  ~TestDenMatGuessCalc() final = default;
  auto calculateGuess() const -> DensityMatrix final {
    DensityMatrix den;
    den.setDensity(Eigen::MatrixXd::Identity(10, 10), 10);
    return den;
  }
};

struct TestInitializer final : public StructureDependentInitializer {
  ~TestInitializer() final = default;

  AtomsOrbitalsIndexes aoIndexes_;
  std::vector<double> coreCharges_;
  void initialize(const Utils::ElementTypeCollection& /*elements*/) {
    aoIndexes_ = AtomsOrbitalsIndexes(10);
    coreCharges_.resize(10);
    for (int i = 0; i < 10; ++i) {
      aoIndexes_.addAtom(1);
      coreCharges_[i] = 1;
    }
  }

  AtomsOrbitalsIndexes getAtomsOrbitalsIndexes() const {
    return aoIndexes_;
  }
  unsigned getNumberElectronsForUnchargedSpecies() const {
    return 10;
  }
  std::vector<double> getCoreCharges() const {
    return coreCharges_;
  }
  bool unrestrictedCalculationPossible() const {
    return false;
  }
};

struct DummyScfMethod : public ScfMethod {
  DummyScfMethod() : ScfMethod(false, Utils::DerivativeOrder::Zero, false) {
    overlapCalculator_ = std::make_shared<TestOverlapCalculator>();
    densityMatrixGuess_ = std::make_shared<TestDenMatGuessCalc>();
    initializer_ = std::make_shared<TestInitializer>();
    electronicPart_ = std::make_shared<TestElectronicPart>(densityMatrix_, aoIndexes_);
    rep_ = std::make_shared<RepulsionCalculator>(elementTypes_, positions_);
    elementTypes_ = ElementTypeCollection(10);
    positions_ = PositionCollection::Zero(10, 3);
    std::fill(elementTypes_.begin(), elementTypes_.end(), ElementType::H);
    for (int i = 0; i < 10; ++i) {
      positions_(i, 0) = i;
    }
    initialize();
  }
  void initialize() final {
    molecularCharge_ = 0;
    ScfMethod::initialize();
  }
};

class LcaoTestMethod : public LcaoMethod {
 public:
  LcaoTestMethod() : LcaoMethod(true, Utils::DerivativeOrder::Two, true) {
  }
  void initialize() final {
    molecularCharge_ = 0;
    nElectronsForUnchargedSpecies_ = 10;
    nAOs_ = 10;
  }
};

struct FieldContribution final : public AdditiveElectronicContribution {
  FieldContribution(const Eigen::MatrixXd& overlap, const Eigen::MatrixXd& mos, const AtomsOrbitalsIndexes& aoIndexes)
    : ovlp_(overlap), mos_(mos), aoIndexes_(aoIndexes) {
  }
  ~FieldContribution() final = default;
  auto isDensityDependent() const -> bool final {
    return true;
  }

  double getElectronicEnergyContribution() const final {
    return energy_;
  }
  // Calculate dipole matrix in MO with mulliken approx and field contribution
  void calculate(const DensityMatrix& densityMatrix, DerivativeOrder /*order*/) final {
    Eigen::MatrixXd positionMatrix(10, 10);
    for (int i = 0; i < 10; ++i) {
      for (int j = 0; j < 10; ++j) {
        positionMatrix(i, j) = i;
      }
    }
    Eigen::MatrixXd popMat = ovlp_.cwiseProduct(densityMatrix.restrictedMatrix());

    // Field in x-direction, 1.0a.u. -> contribution = E * D_x, energy = E * D_x * P
    setElectronicContribution(SpinAdaptedMatrix::createRestricted(0.0001 * (popMat * positionMatrix)));
    energy_ = (popMat * positionMatrix).cwiseProduct(densityMatrix.restrictedMatrix()).sum();
  }

  void addDerivatives(AutomaticDifferentiation::DerivativeContainerType<Derivative::First>& /*derivativeContainer*/) const {
  }
  void addDerivatives(AutomaticDifferentiation::DerivativeContainerType<Derivative::SecondAtomic>& /*derivativeContainer*/) const {
  }
  void addDerivatives(AutomaticDifferentiation::DerivativeContainerType<Derivative::SecondFull>& /*derivativeContainer*/) const {
  }

 private:
  Eigen::MatrixXd ovlp_;
  Eigen::MatrixXd mos_;
  const AtomsOrbitalsIndexes& aoIndexes_;
  double energy_;
};

/**
 * @class ALcaoTest LcaoTest.cpp
 * @brief Comprises tests for the class Scine::Utils::LcaoMethod.
 * @test
 */
class ALcaoTest : public Test {
 public:
  DummyScfMethod scfMethod_;
  LcaoTestMethod lcaoMethod_;

  Core::Log log;

 protected:
  void SetUp() override {
    lcaoMethod_.initialize();
  }
};

TEST_F(ALcaoTest, CorrectlyThrowsInUnrestrictedCalculation) {
  lcaoMethod_.setMolecularCharge(-2);
  lcaoMethod_.setSpinMultiplicity(3);
  lcaoMethod_.setUnrestrictedCalculation(false);

  ASSERT_THROW(lcaoMethod_.verifyPesValidity(), std::runtime_error);
}

TEST_F(ALcaoTest, CorrectlyThrowsInSpinMultiplicity) {
  lcaoMethod_.setMolecularCharge(-2);
  lcaoMethod_.setSpinMultiplicity(2);
  lcaoMethod_.setUnrestrictedCalculation(true);

  ASSERT_THROW(lcaoMethod_.verifyPesValidity(), std::runtime_error);
}

TEST_F(ALcaoTest, DoesNotThrowIfCorrect) {
  lcaoMethod_.setMolecularCharge(-2);
  lcaoMethod_.setSpinMultiplicity(1);
  lcaoMethod_.setUnrestrictedCalculation(true);

  ASSERT_NO_THROW(lcaoMethod_.verifyPesValidity());

  lcaoMethod_.setMolecularCharge(0);
  lcaoMethod_.setSpinMultiplicity(1);
  lcaoMethod_.setUnrestrictedCalculation(false);

  ASSERT_NO_THROW(lcaoMethod_.verifyPesValidity());
}
// TODO: Add more unit tests
TEST_F(ALcaoTest, CanScfAModelSystem) {
  ConvergenceChecker::Thresholds t{1e-9, {}};

  scfMethod_.setConvergenceCriteria(t);
  scfMethod_.setMaxIterations(100);
  Core::Log log;
  scfMethod_.convergedCalculation(log, Utils::Derivative::None);

  EXPECT_LE(*(scfMethod_.getCurrentConvergenceValues()[0]), 1e-9);

  scfMethod_.setConvergenceCriteria({1e-9, 1e-9});
  scfMethod_.setMaxIterations(400);
  scfMethod_.convergedCalculation(log, Utils::Derivative::None);

  EXPECT_LE(*(scfMethod_.getCurrentConvergenceValues()[0]), 1e-9);
  EXPECT_LE(*(scfMethod_.getCurrentConvergenceValues()[1]), 1e-9);
}

TEST_F(ALcaoTest, CanScfAModelSystemWithAdditiveHamiltonian) {
  scfMethod_.setConvergenceCriteria({1e-9, 1e-9});
  scfMethod_.setMaxIterations(800);
  Core::Log log;
  scfMethod_.convergedCalculation(log, Utils::Derivative::None);
  scfMethod_.addElectronicContribution(std::make_shared<FieldContribution>(
      scfMethod_.getOverlapMatrix(), scfMethod_.getMolecularOrbitals().restrictedMatrix(),
      scfMethod_.getAtomsOrbitalsIndexesHolder()));
  scfMethod_.convergedCalculation(log, Utils::Derivative::None);

  EXPECT_LE(*(scfMethod_.getCurrentConvergenceValues()[0]), 1e-9);
  EXPECT_LE(*(scfMethod_.getCurrentConvergenceValues()[1]), 1e-9);
}
} // namespace Tests
} // namespace Utils
} /* namespace Scine */
