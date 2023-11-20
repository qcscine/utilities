/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CrossValidation.h"
#include "Regression/RegressionModel.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <random>

using namespace boost::accumulators;

namespace Scine {
namespace Utils {
namespace MachineLearning {

CrossValidation::CrossValidation(const RegressionModel& model, int k) : model_(model), k_(k) {
  if (k_ <= 1) {
    throw std::runtime_error("The value for k cannot be less than or equal to one.");
  }
}

// Returns MAE and standard deviation as a pair
std::pair<double, double> CrossValidation::evaluateRegressionModel(const Eigen::MatrixXd& featureValues,
                                                                   const Eigen::MatrixXd& targetValues) {
  if (targetValues.rows() != featureValues.rows()) {
    throw std::runtime_error("The number of data points do not match between the feature and target matrices.");
  }

  numDataPoints_ = static_cast<int>(targetValues.rows());
  if (numDataPoints_ % k_ != 0) {
    throw std::runtime_error("The number of data points has to be divisible by k. Set a different value for k.");
  }

  // Shuffle the data randomly
  shuffleData(featureValues, targetValues);

  // Set some important variables
  numDataPointsPerSubset_ = numDataPoints_ / k_;
  numFeatures_ = static_cast<int>(featureValues.cols());
  numTargets_ = static_cast<int>(targetValues.cols());

  std::vector<double> absoluteErrors(k_);

#pragma omp parallel
  {
    auto localModel = model_.clone();
#pragma omp for schedule(dynamic)
    for (int i = 0; i < k_; ++i) {
      performIteration(i, absoluteErrors, *localModel);
    }
  }

  // Calculate MAE and std. dev. and return the result
  return calculateStatistics(absoluteErrors);
}

void CrossValidation::shuffleData(const Eigen::MatrixXd& featureValues, const Eigen::MatrixXd& targetValues) {
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> permutationMatrix(numDataPoints_);
  permutationMatrix.setIdentity();
  std::mt19937 engine(randomSeed_);
  std::shuffle(permutationMatrix.indices().data(),
               permutationMatrix.indices().data() + permutationMatrix.indices().size(), engine);
  shuffledFeatureValues_ = permutationMatrix * featureValues;
  shuffledTargetValues_ = permutationMatrix * targetValues;
}

void CrossValidation::performIteration(int index, std::vector<double>& absoluteErrors, RegressionModel& model) {
  // Create training data set
  // For features:
  Eigen::MatrixXd trainingDataFeatures = transferDataForTraining(shuffledFeatureValues_, index);
  // For targets:
  Eigen::MatrixXd trainingDataTargets = transferDataForTraining(shuffledTargetValues_, index);

  model.trainModel(trainingDataFeatures, trainingDataTargets);

  double errorSum = 0.0;
  // Iterate over the rows that are in the test set
  for (int rowIndex = index * numDataPointsPerSubset_; rowIndex < ((index + 1) * numDataPointsPerSubset_); ++rowIndex) {
    Eigen::VectorXd result = model.predict(shuffledFeatureValues_.row(rowIndex));
    // Iterate over all targets
    for (int targetIndex = 0; targetIndex < numTargets_; ++targetIndex) {
      errorSum += std::abs(result(targetIndex) - shuffledTargetValues_(rowIndex, targetIndex));
    }
  }
  absoluteErrors[index] = errorSum / (numDataPointsPerSubset_ * numTargets_);
}

std::pair<double, double> CrossValidation::calculateStatistics(const std::vector<double>& absoluteErrors) {
  accumulator_set<double, stats<tag::mean, tag::variance(lazy)>> acc;
  std::for_each(absoluteErrors.begin(), absoluteErrors.end(), boost::bind<void>(boost::ref(acc), _1));
  return {mean(acc), sqrt(variance(acc))};
}

void CrossValidation::setNumberOfSubsets(int k) {
  k_ = k;
}

void CrossValidation::setRandomSeed(int seed) {
  randomSeed_ = seed;
}

Eigen::MatrixXd CrossValidation::transferDataForTraining(const Eigen::MatrixXd& originalData, int index) const {
  Eigen::MatrixXd trainingData(shuffledFeatureValues_.rows() - numDataPointsPerSubset_, originalData.cols());
  trainingData.topRows(index * numDataPointsPerSubset_) = originalData.topRows(index * numDataPointsPerSubset_);
  trainingData.bottomRows(trainingData.rows() - index * numDataPointsPerSubset_) =
      originalData.bottomRows(trainingData.rows() - index * numDataPointsPerSubset_);
  return trainingData;
}

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine
