/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_GAUSSIANPROCESSREGRESSION_H
#define UTILS_MATH_GAUSSIANPROCESSREGRESSION_H

#include "../Kernels.h"
#include "RegressionModel.h"
#include <Utils/Constants.h>
#include <Utils/Optimizer/LeastSquares/UpdateFunctionManagerBase.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <boost/optional.hpp>
#include <memory>

namespace Scine {
namespace Utils {
class Settings;
namespace MachineLearning {

/**
 * @class GaussianProcessRegression GaussianProcessRegression.h
 * @brief This class implements Gaussian Processes.
 */
class GaussianProcessRegression : public CloneInterface<GaussianProcessRegression, RegressionModel> {
 public:
  /// @brief Default Constructor.
  GaussianProcessRegression();
  /// @brief Default Destructor.
  ~GaussianProcessRegression() final = default;
  /// @brief Copy Constructor.
  GaussianProcessRegression(const GaussianProcessRegression& rhs);

  struct HyperparameterSpecifier {
    double guess;
    bool toOptimize = true;
    boost::optional<std::pair<double, double>> bounds;
  };

  /**
   * @brief Trains the model applying Gaussian Processes.
   * @param featureValues A matrix containing the input training data which the model is trained on.
   *                      The number of rows is equal to the number of samples and the number of columns is equal
   *                      to the number of features.
   * @param targetValues A matrix containing the target values for the input training data.
   *                     The number of rows is equal to the number of samples and the number of columns is equal
   *                     to the number of targets for each sample (for each data point).
   */
  void trainModel(const Eigen::MatrixXd& featureValues, const Eigen::MatrixXd& targetValues) override;
  /**
   * @brief Predicts the target values for a given data point (set of feature values) with the trained GPR model.
   * @param data A vector containing the values for the features of the data point which shall be predicted by the GPR
   * model.
   * @return A vector of predicted target values. The length of this vector is equal to the number of targets.
   */
  Eigen::VectorXd predict(const Eigen::VectorXd& data) override;
  /**
   * @brief Returns the optimized Hyperparameters.
   */
  Eigen::VectorXd getOptimizedHyperparameters();
  /**
   * @brief Returns the covariance matrix.
   */
  Eigen::MatrixXd getVarianceOfPrediction();

  /**
   * @brief Set the specifications for the hyperparameter theta (guess value, bounds, toOptimize)
   */
  void setThetaHyperparametersGuess(HyperparameterSpecifier specifier);
  /**
   * @brief Set the specifications for the hyperparameter sigma_y^2 (guess value, bounds, toOptimize)
   */
  void setSigmaYSqHyperparametersGuess(HyperparameterSpecifier specifier);
  /**
   * @brief Set the specifications for the hyperparameter sigma_g^2 (guess value, bounds, toOptimize)
   */
  void setSigmaFSqHyperparametersGuess(HyperparameterSpecifier specifier);
  /**
   * @brief Returns the boundaries set for the three hyperparameters.
   */
  auto getBoundaries() const -> std::vector<std::pair<int, std::pair<double, double>>>;

  /**
   * @brief Accessor for the settings.
   * @return Settings& The settings.
   */
  Settings& settings();
  /**
   * @brief Constant accessor for the settings.
   * @return const Settings& The settings.
   */
  const Settings& settings() const;

 private:
  // The kernel function. Default kernel is an autoGaussian kernel. Note that Gaussian processes are currently
  // implemented only for this kernel.
  std::function<AutomaticDifferentiation::FirstND(const Eigen::VectorXd&, const Eigen::VectorXd&, const Kernels::Hyperparameters&)> kernel_ =
      Kernels::autoGaussianKernel;
  /**
   * This functional calculates the negative (log of the) likelihood. Note that it returns the negative of the
   * likelihood value. This is necessary as minimizing the negative likelihood is equal to maximizing the original
   * likelihood value.
   */
  void fit();
  // Constructs the correlation matrix between two sets of data points.
  Eigen::MatrixXd constructCorrelationMatrix(const Eigen::MatrixXd& data1, const Eigen::MatrixXd& data2) const;
  // Target values stored in a matrix (number of target values per data point x number of data points)
  Eigen::MatrixXd targetValues_;
  // Values of the features stored in a matrix (number of features x number of data points)
  Eigen::MatrixXd featureValues_;
  // The number of data points in the training set
  int numDataPoints_ = 0;
  // the number of test points for the prediction
  int numTestPoints_ = 0;
  // Hyperparameters of the kernel
  Kernels::Hyperparameters hyperparameters_;
  // The covariance matrix
  Eigen::MatrixXd variance_;
  // The settings.
  std::unique_ptr<Settings> settings_;
  // The specifier for hyperparameter theta
  HyperparameterSpecifier thetaSpecifier_{
      1.0, true, {{-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()}}};
  // The specifier for hyperparameter sigma_f^2
  HyperparameterSpecifier sigmaFSqSpecifier_{
      1.0, true, {{-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()}}};
  // The specifier for hyperparameter sigma_y^2
  HyperparameterSpecifier sigmaYSqSpecifier_{
      0.1, false, {{-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()}}};
};

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_GAUSSIANPROCESSREGRESSION_H
