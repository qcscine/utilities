/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/Log.h>
#include <Utils/Math/IterativeDiagonalizer/DavidsonDiagonalizer.h>
#include <Utils/Math/IterativeDiagonalizer/DiagonalizerSettings.h>
#include <Utils/Math/IterativeDiagonalizer/IndirectPreconditionerEvaluator.h>
#include <Utils/Math/IterativeDiagonalizer/IndirectSigmaVectorEvaluator.h>
#include <Utils/Math/IterativeDiagonalizer/SubspaceCollapser.h>
#include <Utils/Math/IterativeDiagonalizer/SubspaceOrthogonalizer.h>
#include <gmock/gmock.h>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <chrono>
#include <random>

namespace Scine {
namespace Utils {

using namespace testing;

/**
 * @class ADavidsonDiagonalizerTest @file DavidsonDiagonalizerTest.cpp
 * @test
 */
class ADavidsonDiagonalizerTest : public Test {
 public:
  int arbitraryNumberEigenvalue = 4;
  int arbitraryTrialSubspace = 8;
  int arbitraryDimension = 100;
  Eigen::MatrixXd arbitraryDenseMatrix;
  Eigen::SparseMatrix<double> arbitrarySparseMatrix;
  Core::Log log;

 private:
  void SetUp() final {
    Eigen::VectorXd diagonal = Eigen::VectorXd::Ones(arbitraryDimension);
    for (int i = 0; i < arbitraryDimension; ++i) {
      diagonal(i) *= i;
    }
    Eigen::MatrixXd arbitraryMatrix = 0.01 * Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    arbitraryDenseMatrix = arbitraryMatrix + arbitraryMatrix.transpose();
    arbitraryDenseMatrix.diagonal() += diagonal;

    arbitrarySparseMatrix.resize(arbitraryDimension, arbitraryDimension);

    Eigen::MatrixXd randomMatrix = Eigen::MatrixXd::Random(arbitraryDimension, arbitraryDimension);
    for (Eigen::Index i = 0; i < arbitraryDimension; ++i) {
      for (Eigen::Index j = i; j < arbitraryDimension; ++j) {
        if (i == j || randomMatrix(i, j) > 0.8) { // Random c [-1,1]
          arbitrarySparseMatrix.coeffRef(i, j) = arbitraryDenseMatrix(i, j);
          arbitrarySparseMatrix.coeffRef(j, i) = arbitraryDenseMatrix(i, j);
        }
      }
    }
  }
};

TEST_F(ADavidsonDiagonalizerTest, DavidsonDiagonalizerCanBeConstructed) {
  OrthogonalDavidson davidson{arbitraryNumberEigenvalue, arbitraryDimension};
  davidson.settings().modifyInt(initialGuessDimensionOption, arbitraryTrialSubspace);
  NonOrthogonalDavidson davidsonBalanced{arbitraryNumberEigenvalue, arbitraryDimension};
  davidsonBalanced.settings().modifyInt(initialGuessDimensionOption, arbitraryTrialSubspace);
}

TEST_F(ADavidsonDiagonalizerTest, IndirectPreconditionerEvaluatorCanBeConstructed) {
  Eigen::VectorXd arbitraryDiagonal(5);
  arbitraryDiagonal << 1, 2, 3, 4, 5;
  IndirectPreconditionerEvaluator preconditioner(arbitraryDiagonal);

  // result = -(H - I\lambda)^{-1} * v
  //        = -1/[-5, -4, -3, -2, -1] * v (cwise)
  //        = 1/5, 1/4, 1/3, 1/2, 1
  Eigen::VectorXd expected(5);
  expected << 1. / 5., 1. / 4., 1. / 3., 1. / 2., 1;
  Eigen::VectorXd preconditionedVector = preconditioner.evaluate(Eigen::VectorXd::Ones(5), 6);
  for (int i = 0; i < preconditionedVector.size(); ++i) {
    EXPECT_DOUBLE_EQ(preconditionedVector(i), expected(i));
  }
}

TEST_F(ADavidsonDiagonalizerTest, IndirectSigmaVectorEvaluatorCanBeConstructed) {
  Eigen::MatrixXd arbMatrix(3, 3);
  arbMatrix << 1, 2, 3, 2, 4, 6, 3, 6, 9;
  IndirectSigmaVectorEvaluator<Eigen::MatrixXd> sve(arbMatrix);

  Eigen::MatrixXd guessVectors = Eigen::Vector3d(2, 2, 2).asDiagonal();

  Eigen::MatrixXd result = sve.evaluate(guessVectors);

  ASSERT_EQ(result.rows(), 3);
  ASSERT_EQ(result.cols(), 3);
  for (int row = 0; row < result.rows(); ++row) {
    for (int col = 0; col < result.cols(); ++col) {
      EXPECT_DOUBLE_EQ(result(row, col), arbMatrix(row, col) * 2);
    }
  }

  // collapsed should not have any effect
  sve.collapsed(1);
  result = sve.evaluate(guessVectors);

  ASSERT_EQ(result.rows(), 3);
  ASSERT_EQ(result.cols(), 3);
  for (int row = 0; row < result.rows(); ++row) {
    for (int col = 0; col < result.cols(); ++col) {
      EXPECT_DOUBLE_EQ(result(row, col), arbMatrix(row, col) * 2);
    }
  }

  Eigen::MatrixXd guessVectorsWrongDims = Eigen::Vector2d(2, 2).asDiagonal();
  ASSERT_THROW(sve.evaluate(guessVectorsWrongDims), std::runtime_error);
}

TEST_F(ADavidsonDiagonalizerTest, SparseIndirectSigmaVectorEvaluatorCanBeConstructed) {
  Eigen::MatrixXd arbDenseMatrix(3, 3);
  arbDenseMatrix << 1, 2, 3, 2, 4, 6, 3, 6, 9;
  Eigen::SparseMatrix<double> arbMatrix = arbDenseMatrix.sparseView(0, 0);
  auto sve = std::make_unique<IndirectSigmaVectorEvaluator<Eigen::SparseMatrix<double>>>(arbMatrix);

  Eigen::MatrixXd guessVectors = Eigen::Vector3d(2, 2, 2).asDiagonal();

  Eigen::MatrixXd result = sve->evaluate(guessVectors);

  ASSERT_EQ(result.rows(), 3);
  ASSERT_EQ(result.cols(), 3);
  for (int row = 0; row < result.rows(); ++row) {
    for (int col = 0; col < result.cols(); ++col) {
      EXPECT_DOUBLE_EQ(result(row, col), arbMatrix.coeffRef(row, col) * 2);
    }
  }

  // collapsed should not have any effect
  sve->collapsed(1);
  result = sve->evaluate(guessVectors);

  ASSERT_EQ(result.rows(), 3);
  ASSERT_EQ(result.cols(), 3);
  for (int row = 0; row < result.rows(); ++row) {
    for (int col = 0; col < result.cols(); ++col) {
      EXPECT_DOUBLE_EQ(result(row, col), arbMatrix.coeffRef(row, col) * 2);
    }
  }
  sve.reset();
  ASSERT_FALSE(sve);
}
TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithStandardMethod) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(1);
  auto actualEigenvectors = ref.eigenvectors().leftCols(1);

  OrthogonalDavidson davidson{1, arbitraryDimension};
  davidson.settings().modifyInt(initialGuessDimensionOption, 1);

  davidson.setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(arbitraryDenseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitraryDenseMatrix.diagonal()));
  auto result = davidson.solve(log);

  ASSERT_THAT(result.eigenValues(0), DoubleNear(actualEigenvalues(0), 1e-5));
  for (int j = 0; j < arbitraryDimension; ++j) {
    ASSERT_THAT(abs(result.eigenVectors(j, 0)), DoubleNear(abs(actualEigenvectors(j, 0)), 1e-5));
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithStandardMethodAndManyEigenvaluesWithLargeGuessSpace) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryDimension - 3);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryDimension - 3);

  OrthogonalDavidson davidson{arbitraryDimension - 3, arbitraryDimension};
  davidson.settings().modifyInt(initialGuessDimensionOption, arbitraryDimension - 3);
  davidson.setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(arbitraryDenseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitraryDenseMatrix.diagonal()));

  auto result = davidson.solve(log);

  for (int i = 0; i < arbitraryDimension - 3; ++i) {
    ASSERT_THAT(result.eigenValues(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.eigenVectors(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithStandardMethodAndManyEigenvalues) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryNumberEigenvalue);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryNumberEigenvalue);

  OrthogonalDavidson davidson{arbitraryNumberEigenvalue, arbitraryDimension};
  davidson.settings().modifyInt(initialGuessDimensionOption, arbitraryTrialSubspace);

  davidson.setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(arbitraryDenseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitraryDenseMatrix.diagonal()));

  auto result = davidson.solve(log);

  for (int i = 0; i < arbitraryNumberEigenvalue; ++i) {
    ASSERT_THAT(result.eigenValues(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.eigenVectors(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithSparseMatrix) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitrarySparseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(1);
  auto actualEigenvectors = ref.eigenvectors().leftCols(1);

  OrthogonalDavidson davidson{1, arbitraryDimension};
  davidson.settings().modifyInt(initialGuessDimensionOption, 1);
  davidson.setSigmaVectorEvaluator(
      std::make_unique<IndirectSigmaVectorEvaluator<Eigen::SparseMatrix<double>>>(arbitrarySparseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitrarySparseMatrix.diagonal()));
  auto result = davidson.solve(log);

  ASSERT_THAT(result.eigenValues(0), DoubleNear(actualEigenvalues(0), 1e-5));
  for (int j = 0; j < 100; ++j) {
    ASSERT_THAT(abs(result.eigenVectors(j, 0)), DoubleNear(abs(actualEigenvectors(j, 0)), 1e-5));
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithSparseMatrixAndManyEigenvaluesAndLargeStartingSpace) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitrarySparseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryDimension - 3);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryDimension - 3);

  OrthogonalDavidson davidson{arbitraryDimension - 3, arbitraryDimension};
  davidson.settings().modifyInt(initialGuessDimensionOption, arbitraryDimension - 3);

  davidson.setSigmaVectorEvaluator(
      std::make_unique<IndirectSigmaVectorEvaluator<Eigen::SparseMatrix<double>>>(arbitrarySparseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitrarySparseMatrix.diagonal()));
  auto result = davidson.solve(log);

  for (int i = 0; i < arbitraryDimension - 3; ++i) {
    ASSERT_THAT(result.eigenValues(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.eigenVectors(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithSparseMatrixAndManyEigenvalues) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitrarySparseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryNumberEigenvalue);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryNumberEigenvalue);

  OrthogonalDavidson davidson{arbitraryNumberEigenvalue, arbitraryDimension};
  davidson.settings().modifyInt(initialGuessDimensionOption, arbitraryTrialSubspace);
  davidson.setSigmaVectorEvaluator(
      std::make_unique<IndirectSigmaVectorEvaluator<Eigen::SparseMatrix<double>>>(arbitrarySparseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitrarySparseMatrix.diagonal()));

  auto result = davidson.solve(log);

  for (int i = 0; i < arbitraryNumberEigenvalue; ++i) {
    ASSERT_THAT(result.eigenValues(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryNumberEigenvalue; ++j) {
      ASSERT_THAT(abs(result.eigenVectors(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeMatrixOfSizeOne) {
  Eigen::MatrixXd matrix(1, 1);
  matrix << 13.0;

  OrthogonalDavidson davidson{1, 1};
  davidson.setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(matrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(matrix.diagonal()));

  auto result = davidson.solve(log);

  ASSERT_THAT(result.eigenValues(0), Eq(13.0));
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithSimultDiagBalancedMethod) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(1);
  auto actualEigenvectors = ref.eigenvectors().leftCols(1);

  NonOrthogonalDavidson davidson{1, arbitraryDimension};
  davidson.setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(arbitraryDenseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitraryDenseMatrix.diagonal()));
  davidson.settings().modifyString("gep_algo", "simultaneous_diag");

  auto result = davidson.solve(log);

  EXPECT_THAT(result.eigenValues(0), DoubleNear(actualEigenvalues(0), 1e-5));
  for (int j = 0; j < arbitraryDimension; ++j) {
    EXPECT_THAT(abs(result.eigenVectors(j, 0)), DoubleNear(abs(actualEigenvectors(j, 0)), 1e-5));
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithSimultDiagBalancedMethodAlmostSingular) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(1);
  auto actualEigenvectors = ref.eigenvectors().leftCols(1);
  Eigen::MatrixXd singularGuess(arbitraryDenseMatrix.rows(), 2);
  Eigen::VectorXd guessVector = Eigen::MatrixXd::Identity(arbitraryDenseMatrix.rows(), 1).col(0);
  singularGuess << guessVector + 1e-12 * Eigen::VectorXd::Random(arbitraryDenseMatrix.rows()),
      guessVector + 1e-12 * Eigen::VectorXd::Random(arbitraryDenseMatrix.rows());

  NonOrthogonalDavidson davidson{1, arbitraryDimension};
  davidson.setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(arbitraryDenseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitraryDenseMatrix.diagonal()));
  davidson.settings().modifyInt(initialGuessDimensionOption, 2);
  davidson.setGuess(singularGuess);
  davidson.settings().modifyString("gep_algo", "simultaneous_diag");

  auto result = davidson.solve(log);

  EXPECT_THAT(result.eigenValues(0), DoubleNear(actualEigenvalues(0), 1e-5));
  for (int j = 0; j < arbitraryDimension; ++j) {
    EXPECT_THAT(abs(result.eigenVectors(j, 0)), DoubleNear(abs(actualEigenvectors(j, 0)), 1e-5));
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithStandardBalancedMethod) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(1);
  auto actualEigenvectors = ref.eigenvectors().leftCols(1);

  NonOrthogonalDavidson davidson{1, arbitraryDimension};
  davidson.setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(arbitraryDenseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitraryDenseMatrix.diagonal()));
  davidson.settings().modifyString("gep_algo", "standard");

  auto result = davidson.solve(log);

  EXPECT_THAT(result.eigenValues(0), DoubleNear(actualEigenvalues(0), 1e-5));
  for (int j = 0; j < arbitraryDimension; ++j) {
    EXPECT_THAT(abs(result.eigenVectors(j, 0)), DoubleNear(abs(actualEigenvectors(j, 0)), 1e-5));
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithCholeskyBalancedMethod) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(1);
  auto actualEigenvectors = ref.eigenvectors().leftCols(1);

  NonOrthogonalDavidson davidson{1, arbitraryDimension};
  davidson.setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(arbitraryDenseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitraryDenseMatrix.diagonal()));
  davidson.settings().modifyString("gep_algo", "cholesky");

  auto result = davidson.solve(log);

  EXPECT_THAT(result.eigenValues(0), DoubleNear(actualEigenvalues(0), 1e-5));
  for (int j = 0; j < arbitraryDimension; ++j) {
    EXPECT_THAT(abs(result.eigenVectors(j, 0)), DoubleNear(abs(actualEigenvectors(j, 0)), 1e-5));
  }
}

TEST_F(ADavidsonDiagonalizerTest, ThrowsIfNoCorrectGEPAlgorithmChosen) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);

  NonOrthogonalDavidson davidson{1, arbitraryDimension};
  davidson.setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(arbitraryDenseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitraryDenseMatrix.diagonal()));
  // "vattelappesca" is no valid gep_algo
  davidson.settings().modifyString("gep_algo", "vattelappesca");

  ASSERT_THROW(davidson.solve(log), InvalidDiagonalizerInputException);
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithBalancedStandardMethodAndManyEigenvaluesWithLargeStartingSpace) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryDimension - 3);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryDimension - 3);

  NonOrthogonalDavidson davidson{arbitraryDimension - 3, arbitraryDimension};
  davidson.settings().modifyInt(initialGuessDimensionOption, arbitraryDimension - 3);
  davidson.setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(arbitraryDenseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitraryDenseMatrix.diagonal()));
  auto result = davidson.solve(log);

  for (int i = 0; i < arbitraryDimension - 3; ++i) {
    ASSERT_THAT(result.eigenValues(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.eigenVectors(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithBalancedStandardMethodAndManyEigenvalues) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryNumberEigenvalue);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryNumberEigenvalue);

  NonOrthogonalDavidson davidson{arbitraryNumberEigenvalue, arbitraryDimension};
  davidson.settings().modifyInt(initialGuessDimensionOption, arbitraryTrialSubspace);
  davidson.setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(arbitraryDenseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitraryDenseMatrix.diagonal()));
  auto result = davidson.solve(log);

  for (int i = 0; i < arbitraryNumberEigenvalue; ++i) {
    ASSERT_THAT(result.eigenValues(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.eigenVectors(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithBalancedMethodSparseMatrix) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitrarySparseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(1);
  auto actualEigenvectors = ref.eigenvectors().leftCols(1);

  NonOrthogonalDavidson davidson{1, arbitraryDimension};
  davidson.setSigmaVectorEvaluator(
      std::make_unique<IndirectSigmaVectorEvaluator<Eigen::SparseMatrix<double>>>(arbitrarySparseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitrarySparseMatrix.diagonal()));
  auto result = davidson.solve(log);

  ASSERT_THAT(result.eigenValues(0), DoubleNear(actualEigenvalues(0), 1e-5));
  for (int j = 0; j < 1; ++j) {
    ASSERT_THAT(abs(result.eigenVectors(j, 0)), DoubleNear(abs(actualEigenvectors(j, 0)), 1e-5));
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithBalancedMethodSparseMatrixAndManyEigenvaluesAndLargeStartingSpace) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitrarySparseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryDimension - 3);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryDimension - 3);

  NonOrthogonalDavidson davidson{arbitraryDimension - 3, arbitraryDimension};
  davidson.settings().modifyInt(initialGuessDimensionOption, arbitraryDimension - 3);
  davidson.setSigmaVectorEvaluator(
      std::make_unique<IndirectSigmaVectorEvaluator<Eigen::SparseMatrix<double>>>(arbitrarySparseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitrarySparseMatrix.diagonal()));

  auto result = davidson.solve(log);

  for (int i = 0; i < arbitraryDimension - 3; ++i) {
    ASSERT_THAT(result.eigenValues(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.eigenVectors(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithBalancedMethodSparseMatrixAndManyEigenvalues) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitrarySparseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryNumberEigenvalue);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryNumberEigenvalue);

  NonOrthogonalDavidson davidson{arbitraryNumberEigenvalue, arbitraryDimension};
  davidson.setSigmaVectorEvaluator(
      std::make_unique<IndirectSigmaVectorEvaluator<Eigen::SparseMatrix<double>>>(arbitrarySparseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitrarySparseMatrix.diagonal()));
  auto result = davidson.solve(log);

  for (int i = 0; i < arbitraryNumberEigenvalue; ++i) {
    ASSERT_THAT(result.eigenValues(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.eigenVectors(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, BalancedMethodCanDiagonalizeMatrixOfSizeOne) {
  Eigen::MatrixXd matrix(1, 1);
  matrix << 13.0;

  OrthogonalDavidson davidson{1, 1};
  davidson.setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(matrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(matrix.diagonal()));
  auto result = davidson.solve(log);

  ASSERT_THAT(result.eigenValues(0), Eq(13.0));
}

TEST_F(ADavidsonDiagonalizerTest, SubspaceCollapserWorks) {
  Eigen::VectorXd arbitraryDiagonal(5);

  std::vector<int> notConvergedRoots{0, 1, 2, 3, 4};

  SubspaceCollapser collapser;
  collapser.setMaxSubspaceDimension(10);
  collapser.setEigenvaluesToCompute(5);

  // Subspace eigenvectors for 2 iterations. Only 5 roots are looked for.
  // Dimension of Ritz vectors is 30
  Eigen::MatrixXd eigenvectors1 = Eigen::MatrixXd::Random(30, 5);
  Eigen::MatrixXd eigenvectors2 = Eigen::MatrixXd::Random(30, 5);
  EigenContainer eigenpairs1, eigenpairs2;

  SubspaceOrthogonalizer::qrOrthogonalize(eigenvectors1, 5);
  SubspaceOrthogonalizer::qrOrthogonalize(eigenvectors2, 5);

  eigenpairs1.eigenVectors = eigenvectors1;
  eigenpairs1.eigenValues = arbitraryDiagonal;
  eigenpairs2.eigenVectors = eigenvectors2;
  eigenpairs2.eigenValues = arbitraryDiagonal;

  EXPECT_THAT(collapser.collapseNeeded(eigenpairs1, 9, notConvergedRoots), false);

  EXPECT_THAT(collapser.collapseNeeded(eigenpairs2, 10, notConvergedRoots), true);

  Eigen::MatrixXd collapsed = collapser.getCollapsedOrthogonalSubspace();

  EXPECT_EQ(collapsed.cols(), 10);
  EXPECT_EQ(collapsed.rows(), 30);

  for (int i = 0; i < 30; ++i) {
    for (int j = 0; j < 5; ++j) {
      EXPECT_DOUBLE_EQ(collapsed(i, j), eigenvectors2(i, j));
    }
  }

  // Test that iteration n and n-1 are orthogonal
  for (int i = 0; i < 5; ++i) {
    for (int j = 5; j < 10; ++j) {
      EXPECT_LE(collapsed.col(i).transpose() * collapsed.col(j), 1e-15);
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CorrectlyCalculatesNewIterations) {
  EXPECT_EQ(SubspaceCollapser::calculateSubspaceCollapserIterations(1, 0, 100), 80);
  EXPECT_EQ(SubspaceCollapser::calculateSubspaceCollapserIterations(20, 10, 200), 100);
  EXPECT_EQ(SubspaceCollapser::calculateSubspaceCollapserIterations(100, 0, 1000), 500);
  EXPECT_EQ(SubspaceCollapser::calculateSubspaceCollapserIterations(100, 0, 100), 100);
}

TEST_F(ADavidsonDiagonalizerTest, NonOrthogonalSubspaceCollapserWorks) {
  Eigen::VectorXd arbitraryDiagonal(5);

  std::vector<int> notConvergedRoots{0, 1, 2, 3, 4};

  SubspaceCollapser collapser;
  collapser.setMaxSubspaceDimension(10);
  collapser.setEigenvaluesToCompute(5);

  // Subspace eigenvectors for 2 iterations. Only 5 roots are looked for.
  // Dimension of Ritz vectors is 30
  Eigen::MatrixXd eigenvectors1 = Eigen::MatrixXd::Random(30, 5);
  Eigen::MatrixXd eigenvectors2 = Eigen::MatrixXd::Random(30, 5);
  EigenContainer eigenpairs1, eigenpairs2;

  SubspaceOrthogonalizer::qrOrthogonalize(eigenvectors1, 5);
  SubspaceOrthogonalizer::qrOrthogonalize(eigenvectors2, 5);

  eigenpairs1.eigenVectors = eigenvectors1;
  eigenpairs1.eigenValues = arbitraryDiagonal;
  eigenpairs2.eigenVectors = eigenvectors2;
  eigenpairs2.eigenValues = arbitraryDiagonal;

  EXPECT_THAT(collapser.collapseNeeded(eigenpairs1, 9, notConvergedRoots), false);

  EXPECT_THAT(collapser.collapseNeeded(eigenpairs2, 10, notConvergedRoots), true);

  Eigen::MatrixXd collapsed = collapser.getCollapsedNonOrthogonalSubspace();

  EXPECT_EQ(collapsed.cols(), 10);
  EXPECT_EQ(collapsed.rows(), 30);

  for (int i = 0; i < 30; ++i) {
    for (int j = 0; j < 5; ++j) {
      EXPECT_THAT(collapsed(i, j), DoubleNear(eigenvectors2(i, j), 1e-10));
    }
  }

  // Test that iteration n and n-1 are orthogonal
  for (int i = 0; i < 5; ++i) {
    for (int j = 5; j < 10; ++j) {
      EXPECT_LE(collapsed.col(i).transpose() * collapsed.col(j), 1e-15);
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDoDavidsonWithSubsetCollapse) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryNumberEigenvalue);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryNumberEigenvalue);

  OrthogonalDavidson davidson{arbitraryNumberEigenvalue, arbitraryDimension};
  davidson.settings().modifyInt(initialGuessDimensionOption, arbitraryTrialSubspace);
  davidson.setSigmaVectorEvaluator(std::make_unique<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(arbitraryDenseMatrix));
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitraryDenseMatrix.diagonal()));
  davidson.settings().modifyInt(subspaceCollapseDimensionOption, 10);

  auto result = davidson.solve(log);

  for (int i = 0; i < arbitraryNumberEigenvalue; ++i) {
    ASSERT_THAT(result.eigenValues(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.eigenVectors(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
  auto result2 = davidson.getEigenPairs();
  for (int i = 0; i < arbitraryNumberEigenvalue; ++i) {
    ASSERT_THAT(result.eigenValues(i), DoubleNear(result2.eigenValues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.eigenVectors(j, i)), DoubleNear(abs(result2.eigenVectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, IterativeDiagonalizerUtilitiesFunctions) {
  OrthogonalDavidson davidson{arbitraryNumberEigenvalue, arbitraryDimension};
  davidson.settings().modifyInt(initialGuessDimensionOption, arbitraryTrialSubspace);
  davidson.setPreconditionerEvaluator(std::make_unique<IndirectPreconditionerEvaluator>(arbitraryDenseMatrix.diagonal()));
  davidson.settings().modifyInt(subspaceCollapseDimensionOption, 10);
  ASSERT_THROW(davidson.solve(log), InvalidSigmaVectorEvaluator);

  auto sve = std::make_shared<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(arbitraryDenseMatrix);
  davidson.setSigmaVectorEvaluator(sve);

  auto& sve2 = davidson.getSigmaVectorEvaluator();
  ASSERT_EQ(sve.get(), &sve2);

  davidson.settings().modifyInt(SettingsNames::maxDavidsonIterations, 0);
  ASSERT_THROW(davidson.solve(log), DiagonalizerNotConvergedException);

  InvalidSigmaVectorEvaluator testException;
  DiagonalizerNotConvergedException testException2;

  ASSERT_EQ(testException.what(),
            "Empty sigma vector evaluator or preconditioner for Davidson iterative diagonalizer.");
  ASSERT_EQ(testException2.what(), "Davidson could not converge.");

  const auto& setting = davidson.settings();
  ASSERT_EQ(setting.getInt(subspaceCollapseDimensionOption), davidson.settings().getInt(subspaceCollapseDimensionOption));

  ASSERT_THROW(davidson.setGuess(Eigen::MatrixXd(Eigen::MatrixXd::Identity(200, 200))), std::runtime_error);
}

TEST_F(ADavidsonDiagonalizerTest, TestDavidsonSettingsFunctions) {
  auto testException = InvalidDiagonalizerInputException("This is a test exception");
  ASSERT_EQ(std::string(testException.what()), "Input error: This is a test exception");
  auto settings = DiagonalizerSettings(100, 1000);
  settings.modifyInt(numberOfRootsOption, -1);
  ASSERT_THROW(settings.check(1000), InvalidDiagonalizerInputException);
  ASSERT_THROW(DiagonalizerSettings(1000, 100), InvalidDiagonalizerInputException);
  auto settings2 = DiagonalizerSettings(100, 1000);
  settings2.modifyInt(initialGuessDimensionOption, 2000);
  ASSERT_THROW(settings2.check(1000), InvalidDiagonalizerInputException);
}

} // namespace Utils
} // namespace Scine
