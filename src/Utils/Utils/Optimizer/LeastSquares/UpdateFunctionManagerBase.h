/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_UPDATEFUNCTIONMANAGERBASE_H_
#define UTILS_UPDATEFUNCTIONMANAGERBASE_H_

#include <Eigen/Core>

namespace Scine {
namespace Utils {
/**
 * @class UpdateFunctionManagerBase UpdateFunctionManagerBase.h
 * @brief This class manages the two different update functions that are neccessary for
 *        a least squares optimization using Eigen.
 */
class UpdateFunctionManagerBase {
 public:
  /**
   * @brief Update the errors vector for the least squares optimization for a given set of parameters.
   * @param parameters The parameters.
   * @param errors The error vector to update within this function.
   */
  virtual void updateErrors(const Eigen::VectorXd& parameters, Eigen::VectorXd& errors) = 0;
  /**
   * @brief Update the Jacobian matrix for the least squares optimization for a given set of parameters.
   *        The base implementation is a numerical Jacobian.
   * @param parameters The parameters.
   * @param jacobian The Jacobian matrix to update within this function.
   */
  virtual void updateJacobian(const Eigen::VectorXd& parameters, Eigen::MatrixXd& jacobian);
  /**
   * @brief This function returns the number of data points present in the least squares optimization.
   * @param parameters The parameters.
   * @return The number of data points present in the least squares optimization.
   */
  virtual int getNumberOfDataPoints(const Eigen::VectorXd& parameters) const = 0;
};

// Base implementation of a numerical Jacobian
inline void UpdateFunctionManagerBase::updateJacobian(const Eigen::VectorXd& parameters, Eigen::MatrixXd& jacobian) {
  constexpr double oneThird = 1.0 / 3.0;
  int numberOfDataPoints = getNumberOfDataPoints(parameters);
  auto numberOfParameters = static_cast<int>(parameters.size());
  jacobian.resize(numberOfDataPoints, numberOfParameters);
  jacobian.setZero();

  double stepsize;

  const double epsilon = std::numeric_limits<double>::epsilon(); // optimal step size

  Eigen::VectorXd copyOfParameters = parameters; // make a copy

  for (int i = 0; i < numberOfParameters; ++i) {
    auto x_temp = copyOfParameters(i);

    if (std::abs(x_temp) < 0.05) { // small values of x
      stepsize = 1e-5;
    }
    else {
      stepsize = x_temp * pow(epsilon, oneThird); // other values of x
    }

    copyOfParameters(i) = x_temp + stepsize;
    Eigen::VectorXd errorsPlus;
    updateErrors(copyOfParameters, errorsPlus);

    copyOfParameters(i) = x_temp - stepsize;
    Eigen::VectorXd errorsMinus;
    updateErrors(copyOfParameters, errorsMinus);

    copyOfParameters(i) = x_temp;
    Eigen::VectorXd grad = (errorsPlus - errorsMinus) / (2 * stepsize);

    // against numerical instabilities
    for (int j = 0; j < grad.size(); ++j) {
      if (std::abs(grad(j)) < 1e-8)
        grad(j) = 0;
    }

    jacobian.block(0, i, numberOfDataPoints, 1) = grad;
  }
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_UPDATEFUNCTIONMANAGERBASE_H_
