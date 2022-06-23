/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_REGRESSIONMODEL_H
#define UTILS_MATH_REGRESSIONMODEL_H

#include <Eigen/Dense>
#include <memory>

namespace Scine {
namespace Utils {
class Settings;
namespace MachineLearning {

/**
 * @class RegressionModel RegressionModel.h
 * @brief Base class for a regression model.
 */
class RegressionModel {
 public:
  /// @brief Default constructor.
  RegressionModel() = default;
  /// @brief Default destructor.
  virtual ~RegressionModel() = default;
  /**
   * @brief Trains the model.
   * @param featureValues A matrix containing the input training data which the model is trained on.
   *                      The number of rows is equal to the number of samples and the number of columns is equal
   *                      to the number of features.
   * @param targetValues A matrix containing the target values for the input training data.
   *                     The number of rows is equal to the number of samples and the number of columns is equal
   *                     to the number of targets for each sample (for each data point).
   */
  virtual void trainModel(const Eigen::MatrixXd& featureValues, const Eigen::MatrixXd& targetValues) = 0;
  /**
   * @brief Predicts the target values for a given data point (set of feature values).
   * @param data A vector containing the values for the features of the data point which shall be predicted by the model.
   * @return A vector of predicted target values. The length of this vector is equal to the number of targets.
   */
  virtual Eigen::VectorXd predict(const Eigen::VectorXd& data) = 0;

  /**
   * @brief Method allowing to clone the derived class into a RegressionModel
   *        The derived, leaf class needs to inherit from Utils::CloneInterface and,
   *        if needed, implement a custom copy constructor. This reduces boilerplate code.
   */
  std::unique_ptr<RegressionModel> clone() const {
    return std::unique_ptr<RegressionModel>(this->cloneImpl());
  }

 private:
  /*
   * Implementation of the clone() function, pure virtual private method.
   * It returns a pointer to allow for covariant return types in inheritance.
   */
  virtual RegressionModel* cloneImpl() const = 0;
};

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_REGRESSIONMODEL_H
