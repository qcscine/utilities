/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_CROSSVALIDATION_H
#define UTILS_MATH_CROSSVALIDATION_H

#include <Eigen/Dense>
#include <vector>

namespace Scine {
namespace Utils {
namespace MachineLearning {
class RegressionModel;

/**
 * @class CrossValidation CrossValidation.h
 * @brief This class allows for the estimation of the generalization error of a model via k-fold cross validation.
 */
class CrossValidation {
 public:
  /**
   * @brief Constructor.
   * @param model The regression model.
   * @param k The number of subsets k. If k is equal to the number of data points, the evaluation
   *          becomes a leave-one-out cross validation.
   */
  CrossValidation(const RegressionModel& model, int k);
  /**
   * @brief Calculates the k-fold cross validation of a given data set with a given regression model.
   * @param featureValues The feature matrix. The number of rows is equal to the number of samples. Columns = features.
   * @param targetValues The target matrix. The number of rows is equal to the number of samples. Columns = targets.
   * @return The mean absolute error and the standard deviation of the model evaluation as a std::pair.
   */
  std::pair<double, double> evaluateRegressionModel(const Eigen::MatrixXd& featureValues, const Eigen::MatrixXd& targetValues);
  /**
   * @brief Setter for the number of folds k.
   */
  void setNumberOfSubsets(int k);
  /**
   * @brief Setter for the random seed.
   */
  void setRandomSeed(int seed);

 private:
  /*
   * @brief Shuffles the data randomly.
   */
  void shuffleData(const Eigen::MatrixXd& featureValues, const Eigen::MatrixXd& targetValues);
  /*
   * @brief Performs one iteration (with the index 'index') of the cross validation algorithm (train model and evaluate)
   */
  void performIteration(int index, std::vector<double>& absoluteErrors, RegressionModel& model);
  /*
   * @brief Transfers the data from the original feature and target matrix to one for training and returns it.
   */
  Eigen::MatrixXd transferDataForTraining(const Eigen::MatrixXd& originalData, int index) const;
  /*
   * @brief Calculates the mean and the standard deviation of the given vector of absolute errors
   */
  static std::pair<double, double> calculateStatistics(const std::vector<double>& absoluteErrors);
  // Regression model
  const RegressionModel& model_;
  // Number of subsets k.
  int k_;
  // Random seed.
  int randomSeed_ = 42;
  // The shuffled feature matrix
  Eigen::MatrixXd shuffledFeatureValues_;
  // The shuffled target matrix
  Eigen::MatrixXd shuffledTargetValues_;
  // Number of data points
  int numDataPoints_ = 0;
  // Number of data points in each subset
  int numDataPointsPerSubset_ = 0;
  // Number of features of the data
  int numFeatures_ = 0;
  // Number of targets in the data set
  int numTargets_ = 0;
};

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_CROSSVALIDATION_H
