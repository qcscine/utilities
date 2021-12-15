/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/IO/FormattedIOUtils.h>
#include <gmock/gmock.h>
#include <boost/filesystem.hpp>
#include <fstream>

namespace Scine {
namespace Utils {
using namespace testing;

/**
 * @class FormattedIOUtilsTest FormattedIOUtilsTest.cpp
 * @brief Tests the io functionalities for converting Eigen matrices to and from csv files.
 * @test
 */
class FormattedIOUtilsTest : public Test {
 public:
  void TearDown() override {
    boost::filesystem::remove_all(testFilename);
  }
  const char* testFilename = "test_formatted_io_utils.csv";
};

TEST_F(FormattedIOUtilsTest, MatrixToAndFromCsvWithCommaDelimiter) {
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Random(32, 5);

  std::ofstream file(testFilename);
  matrixToCsv(file, matrix, ',');
  file.close();
  ASSERT_TRUE(boost::filesystem::exists(testFilename));

  Eigen::MatrixXd readMatrix = csvToMatrix(testFilename, ',');
  ASSERT_THAT(readMatrix.size(), Eq(matrix.size()));
  ASSERT_THAT(readMatrix.rows(), Eq(matrix.rows()));
  ASSERT_THAT(readMatrix.cols(), Eq(matrix.cols()));
  ASSERT_TRUE(readMatrix.isApprox(matrix, 1e-6));
}

TEST_F(FormattedIOUtilsTest, MatrixToAndFromCsvWithSpaceDelimiter) {
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Random(4, 4);

  std::ofstream file(testFilename);
  matrixToCsv(file, matrix, ' ');
  file.close();
  ASSERT_TRUE(boost::filesystem::exists(testFilename));

  Eigen::MatrixXd readMatrix = csvToMatrix(testFilename, ' ');
  ASSERT_THAT(readMatrix.size(), Eq(matrix.size()));
  ASSERT_THAT(readMatrix.rows(), Eq(matrix.rows()));
  ASSERT_THAT(readMatrix.cols(), Eq(matrix.cols()));
  ASSERT_TRUE(readMatrix.isApprox(matrix, 1e-6));
}

TEST_F(FormattedIOUtilsTest, ColumnVectorToAndFromCsv) {
  Eigen::VectorXd vector = Eigen::VectorXd::Random(12);

  std::ofstream file(testFilename);
  matrixToCsv(file, vector, ':');
  file.close();
  ASSERT_TRUE(boost::filesystem::exists(testFilename));

  Eigen::MatrixXd readMatrix = csvToMatrix(testFilename, ':');
  ASSERT_THAT(readMatrix.size(), Eq(vector.size()));
  ASSERT_THAT(readMatrix.rows(), Eq(vector.size()));
  ASSERT_THAT(readMatrix.cols(), Eq(1));
  ASSERT_TRUE(readMatrix.isApprox(vector, 1e-6));
}

TEST_F(FormattedIOUtilsTest, RowVectorToAndFromCsv) {
  Eigen::RowVectorXd vector = Eigen::RowVectorXd::Random(12);

  std::ofstream file(testFilename);
  matrixToCsv(file, vector, ':');
  file.close();
  ASSERT_TRUE(boost::filesystem::exists(testFilename));

  Eigen::MatrixXd readMatrix = csvToMatrix(testFilename, ':');
  ASSERT_THAT(readMatrix.size(), Eq(vector.size()));
  ASSERT_THAT(readMatrix.rows(), Eq(1));
  ASSERT_THAT(readMatrix.cols(), Eq(vector.size()));
  ASSERT_TRUE(readMatrix.isApprox(vector, 1e-6));
}

} // namespace Utils
} // namespace Scine
