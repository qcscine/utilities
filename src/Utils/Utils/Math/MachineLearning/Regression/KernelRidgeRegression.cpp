/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "KernelRidgeRegression.h"
#include "RegressionSettings.h"

namespace Scine {
namespace Utils {
namespace MachineLearning {

void KernelRidgeRegression::trainModel(const Eigen::MatrixXd& featureValues, const Eigen::MatrixXd& targetValues) {
  if (targetValues.rows() != featureValues.rows()) {
    throw std::runtime_error("The number of data points do not match between the feature and target matrices.");
  }

  // Transpose the feature and target matrices since they are column-major.
  targetValues_ = targetValues.transpose();
  featureValues_ = featureValues.transpose();
  numDataPoints_ = static_cast<int>(targetValues_.cols());

  // Construct the kernel matrix
  Eigen::MatrixXd kernelMatrix(numDataPoints_, numDataPoints_);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < numDataPoints_; ++i) {
    for (int j = i; j < numDataPoints_; ++j) {
      double value = kernel_(featureValues_.col(i), featureValues_.col(j), hyperparameters_);
      kernelMatrix(i, j) = value;
    }
  }

  kernelMatrix = kernelMatrix.selfadjointView<Eigen::Upper>();

  // Inversion of the regularized kernel matrix
  Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(kernelMatrix.rows(), kernelMatrix.cols());
  Eigen::MatrixXd regularizedKernelMatrix = kernelMatrix + regularizationFactor_ * identity;
  invertedRegularizedKernelMatrix_ = regularizedKernelMatrix.inverse();
}

Eigen::VectorXd KernelRidgeRegression::predict(const Eigen::VectorXd& data) {
  if (invertedRegularizedKernelMatrix_.size() == 0) {
    throw std::runtime_error("The model has not been trained yet!");
  }

  // Calculate kernel vector
  Eigen::VectorXd kernelVector(numDataPoints_);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < numDataPoints_; ++i) {
    kernelVector(i) = kernel_(featureValues_.col(i), data, hyperparameters_);
  }

  // Calculate and return prediction
  Eigen::MatrixXd kernelProduct = invertedRegularizedKernelMatrix_ * kernelVector;
  Eigen::VectorXd prediction = targetValues_ * kernelProduct;

  for (int k = 0; k < prediction.size(); ++k) {
    if (std::isnan(prediction(k))) {
      throw std::runtime_error(
          "One of the predicted results could not be obtained. The reason is probably that the regularization factor "
          "is too small and therefore causes problems during the inversion of the kernel matrix.");
    }
  }

  return prediction;
}

void KernelRidgeRegression::setRegularizationFactor(double factor) {
  regularizationFactor_ = factor;
}

void KernelRidgeRegression::setKernel(std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&)> kernel,
                                      Eigen::VectorXd hyperparameters) {
  kernel_ = std::move(kernel);
  hyperparameters_ = std::move(hyperparameters);
}

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine
