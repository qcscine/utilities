/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/IO/Logger.h>
#include <Utils/Math/IterativeDiagonalizer/DavidsonDiagonalizer.h>
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
        double randomValue = arbitraryDenseMatrix(i, j);
        if (i == j)
          arbitrarySparseMatrix.coeffRef(i, j) = arbitraryDenseMatrix(i, j);
        else if (randomMatrix(i, j) > 0.8) // Random c [-1,1]
          arbitrarySparseMatrix.coeffRef(i, j) = arbitraryDenseMatrix(i, j);
      }
    }
  }
};

TEST_F(ADavidsonDiagonalizerTest, DavidsonDiagonalizerCanBeConstructed) {
  DavidsonDiagonalizer<Eigen::MatrixXd> davidson{arbitraryNumberEigenvalue, arbitraryTrialSubspace, arbitraryDimension};
  DavidsonDiagonalizer<Eigen::SparseMatrix<double>> davidsonSparse{arbitraryNumberEigenvalue, arbitraryTrialSubspace,
                                                                   arbitraryDimension};
  DavidsonDiagonalizer<Eigen::MatrixXd, DavidsonBalancedType::balanced> davidsonBalanced{
      arbitraryNumberEigenvalue, arbitraryTrialSubspace, arbitraryDimension};
  DavidsonDiagonalizer<Eigen::SparseMatrix<double>, DavidsonBalancedType::balanced> davidsonSparseBalanced{
      arbitraryNumberEigenvalue, arbitraryTrialSubspace, arbitraryDimension};
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithStandardMethod) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(1);
  auto actualEigenvectors = ref.eigenvectors().leftCols(1);

  DavidsonDiagonalizer<Eigen::MatrixXd> davidson{1, 1, arbitraryDimension};
  davidson.setMatrixToDiagonalize(arbitraryDenseMatrix);
  auto result = davidson.solve();

  ASSERT_THAT(result.first(0), DoubleNear(actualEigenvalues(0), 1e-5));
  for (int j = 0; j < arbitraryDimension; ++j) {
    ASSERT_THAT(abs(result.second(j, 0)), DoubleNear(abs(actualEigenvectors(j, 0)), 1e-5));
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithStandardMethodAndManyEigenvaluesWithLargeGuessSpace) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryDimension - 3);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryDimension - 3);

  DavidsonDiagonalizer<Eigen::MatrixXd> davidson{arbitraryDimension - 3, arbitraryDimension - 3, arbitraryDimension};
  davidson.setMatrixToDiagonalize(arbitraryDenseMatrix);
  auto result = davidson.solve();

  for (int i = 0; i < arbitraryDimension - 3; ++i) {
    ASSERT_THAT(result.first(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.second(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithStandardMethodAndManyEigenvalues) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryNumberEigenvalue);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryNumberEigenvalue);

  DavidsonDiagonalizer<Eigen::MatrixXd> davidson{arbitraryNumberEigenvalue, arbitraryTrialSubspace, arbitraryDimension};
  davidson.setMatrixToDiagonalize(arbitraryDenseMatrix);
  auto result = davidson.solve();

  for (int i = 0; i < arbitraryNumberEigenvalue; ++i) {
    ASSERT_THAT(result.first(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.second(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithSparseMatrix) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitrarySparseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(1);
  auto actualEigenvectors = ref.eigenvectors().leftCols(1);

  DavidsonDiagonalizer<Eigen::SparseMatrix<double>> davidson{1, 1, arbitraryDimension};
  davidson.setMatrixToDiagonalize(arbitrarySparseMatrix);
  auto result = davidson.solve();

  ASSERT_THAT(result.first(0), DoubleNear(actualEigenvalues(0), 1e-5));
  for (int j = 0; j < 1; ++j) {
    ASSERT_THAT(abs(result.second(j, 0)), DoubleNear(abs(actualEigenvectors(j, 0)), 1e-5));
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithSparseMatrixAndManyEigenvaluesAndLargeStartingSpace) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitrarySparseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryDimension - 3);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryDimension - 3);

  DavidsonDiagonalizer<Eigen::SparseMatrix<double>> davidson{arbitraryDimension - 3, arbitraryDimension - 3, arbitraryDimension};
  davidson.setMatrixToDiagonalize(arbitrarySparseMatrix);
  auto result = davidson.solve();

  for (int i = 0; i < arbitraryDimension - 3; ++i) {
    ASSERT_THAT(result.first(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.second(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithSparseMatrixAndManyEigenvalues) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitrarySparseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryNumberEigenvalue);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryNumberEigenvalue);

  DavidsonDiagonalizer<Eigen::SparseMatrix<double>> davidson{arbitraryNumberEigenvalue, arbitraryTrialSubspace,
                                                             arbitraryDimension};
  davidson.setMatrixToDiagonalize(arbitrarySparseMatrix);
  auto result = davidson.solve();

  for (int i = 0; i < arbitraryNumberEigenvalue; ++i) {
    ASSERT_THAT(result.first(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryNumberEigenvalue; ++j) {
      ASSERT_THAT(abs(result.second(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeMatrixOfSizeOne) {
  Eigen::MatrixXd matrix(1, 1);
  matrix << 13.0;

  DavidsonDiagonalizer<Eigen::MatrixXd> davidson{1, 1, 1};
  davidson.setMatrixToDiagonalize(matrix);
  auto result = davidson.solve();

  ASSERT_THAT(result.first(0), Eq(13.0));
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithStandardBalancedMethod) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(1);
  auto actualEigenvectors = ref.eigenvectors().leftCols(1);

  DavidsonDiagonalizer<Eigen::MatrixXd, DavidsonBalancedType::balanced> davidson{1, 1, arbitraryDimension};
  davidson.setMatrixToDiagonalize(arbitraryDenseMatrix);
  auto result = davidson.solve();

  EXPECT_THAT(result.first(0), DoubleNear(actualEigenvalues(0), 1e-5));
  for (int j = 0; j < arbitraryDimension; ++j) {
    EXPECT_THAT(abs(result.second(j, 0)), DoubleNear(abs(actualEigenvectors(j, 0)), 1e-5));
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithBalancedStandardMethodAndManyEigenvaluesWithLargeStartingSpace) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryDimension - 3);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryDimension - 3);

  DavidsonDiagonalizer<Eigen::MatrixXd, DavidsonBalancedType::balanced> davidson{
      arbitraryDimension - 3, arbitraryDimension - 3, arbitraryDimension};
  davidson.setMatrixToDiagonalize(arbitraryDenseMatrix);
  auto result = davidson.solve();

  for (int i = 0; i < arbitraryDimension - 3; ++i) {
    ASSERT_THAT(result.first(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.second(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithBalancedStandardMethodAndManyEigenvalues) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitraryDenseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryNumberEigenvalue);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryNumberEigenvalue);

  DavidsonDiagonalizer<Eigen::MatrixXd, DavidsonBalancedType::balanced> davidson{
      arbitraryNumberEigenvalue, arbitraryTrialSubspace, arbitraryDimension};
  davidson.setMatrixToDiagonalize(arbitraryDenseMatrix);
  auto result = davidson.solve();

  for (int i = 0; i < arbitraryNumberEigenvalue; ++i) {
    ASSERT_THAT(result.first(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.second(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithBalancedMethodSparseMatrix) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitrarySparseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(1);
  auto actualEigenvectors = ref.eigenvectors().leftCols(1);

  DavidsonDiagonalizer<Eigen::SparseMatrix<double>, DavidsonBalancedType::balanced> davidson{1, 1, arbitraryDimension};
  davidson.setMatrixToDiagonalize(arbitrarySparseMatrix);
  auto result = davidson.solve();

  ASSERT_THAT(result.first(0), DoubleNear(actualEigenvalues(0), 1e-5));
  for (int j = 0; j < 1; ++j) {
    ASSERT_THAT(abs(result.second(j, 0)), DoubleNear(abs(actualEigenvectors(j, 0)), 1e-5));
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithBalancedMethodSparseMatrixAndManyEigenvaluesAndLargeStartingSpace) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitrarySparseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryDimension - 3);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryDimension - 3);

  DavidsonDiagonalizer<Eigen::SparseMatrix<double>, DavidsonBalancedType::balanced> davidson{
      arbitraryDimension - 3, arbitraryDimension - 3, arbitraryDimension};
  davidson.setMatrixToDiagonalize(arbitrarySparseMatrix);
  auto result = davidson.solve();

  for (int i = 0; i < arbitraryDimension - 3; ++i) {
    ASSERT_THAT(result.first(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.second(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, CanDiagonalizeWithBalancedMethodSparseMatrixAndManyEigenvalues) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(arbitrarySparseMatrix);
  auto actualEigenvalues = ref.eigenvalues().head(arbitraryNumberEigenvalue);
  auto actualEigenvectors = ref.eigenvectors().leftCols(arbitraryNumberEigenvalue);

  DavidsonDiagonalizer<Eigen::SparseMatrix<double>, DavidsonBalancedType::balanced> davidson{
      arbitraryNumberEigenvalue, arbitraryNumberEigenvalue, arbitraryDimension};
  davidson.setMatrixToDiagonalize(arbitrarySparseMatrix);
  auto result = davidson.solve();

  for (int i = 0; i < arbitraryNumberEigenvalue; ++i) {
    ASSERT_THAT(result.first(i), DoubleNear(actualEigenvalues(i), 1e-5));
    for (int j = 0; j < arbitraryDimension; ++j) {
      ASSERT_THAT(abs(result.second(j, i)), DoubleNear(abs(actualEigenvectors(j, i)), 1e-5));
    }
  }
}

TEST_F(ADavidsonDiagonalizerTest, BalancedMethodCanDiagonalizeMatrixOfSizeOne) {
  Eigen::MatrixXd matrix(1, 1);
  matrix << 13.0;

  DavidsonDiagonalizer<Eigen::MatrixXd> davidson{1, 1, 1};
  davidson.setMatrixToDiagonalize(matrix);
  auto result = davidson.solve();

  ASSERT_THAT(result.first(0), Eq(13.0));
}
} // namespace Utils
} // namespace Scine
