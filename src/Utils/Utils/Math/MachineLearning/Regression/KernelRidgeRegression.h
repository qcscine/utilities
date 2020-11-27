/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_KERNELRIDGEREGRESSION_H
#define UTILS_MATH_KERNELRIDGEREGRESSION_H

#include "../Kernels.h"
#include "RegressionModel.h"
#include <Utils/Technical/CloneInterface.h>
#include <memory>

namespace Scine {
namespace Utils {
namespace MachineLearning {

/**
 * @class KernelRidgeRegression KernelRidgeRegression.h
 * @brief This class implements the kernelized version of the Ridge Regression method.
 */
class KernelRidgeRegression : public CloneInterface<KernelRidgeRegression, RegressionModel> {
 public:
  /**
   * @brief Trains the model applying Kernel Ridge Regression.
   * @param featureValues A matrix containing the input training data which the model is trained on.
   *                      The number of rows is equal to the number of samples and the number of columns is equal
   *                      to the number of features.
   * @param targetValues A matrix containing the target values for the input training data.
   *                     The number of rows is equal to the number of samples and the number of columns is equal
   *                     to the number of targets for each sample (for each data point).
   */
  void trainModel(const Eigen::MatrixXd& featureValues, const Eigen::MatrixXd& targetValues) override;
  /**
   * @brief Predicts the target values for a given data point (set of feature values) with the trained KRR model.
   * @param data A vector containing the values for the features of the data point which shall be predicted by the KRR
   * model.
   * @return A vector of predicted target values. The length of this vector is equal to the number of targets.
   */
  Eigen::VectorXd predict(const Eigen::VectorXd& data) const override;

  /**
   * @brief Setter for the regularization factor. Note that a very small value can lead to overfitting or to
   *        problems during the inversion of the kernel matrix. It has to be set before training the KRR model.
   * @param factor The regularization factor.
   */
  void setRegularizationFactor(double factor);
  /**
   * @brief Setter for the kernel and its hyperparameters.
   * @param kernel The kernel function.
   * @param hyperparameters The hyperparameters of the kernel function that was provided. For the meaning of the
   *                        hyperparameters of the standard kernels, see the Kernels.h file, where some standard
   *                        kernels are implemented.
   */
  void setKernel(std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&, const std::vector<double>&)> kernel,
                 std::vector<double> hyperparameters);

 private:
  // The kernel function.
  std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&, const std::vector<double>&)> kernel_ =
      Kernels::linearKernel;
  // The inverted regularized kernel matrix: (K + alpha*I)^(-1)
  Eigen::MatrixXd invertedRegularizedKernelMatrix_;
  // Target values stored in a matrix (number of target values per data point x number of data points)
  Eigen::MatrixXd targetValues_;
  // Values of the features stored in a matrix (number of features x number of data points)
  Eigen::MatrixXd featureValues_;
  // Number of data points
  int numDataPoints_ = 0;
  // Regularization factor
  double regularizationFactor_ = 1e-3;
  // Hyperparameters of the kernel
  std::vector<double> hyperparameters_;
};

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_KERNELRIDGEREGRESSION_H
