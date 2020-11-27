/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LevenbergMarquardt.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>

namespace Scine {
namespace Utils {

void LevenbergMarquardt::optimize(Eigen::VectorXd& parameters, UpdateFunctionManagerBase& updateFunctionManager) {
  LMFunctor functor(updateFunctionManager);
  functor.n = static_cast<int>(parameters.size());
  functor.m = functor.updateFunctionManager_.getNumberOfDataPoints(parameters);
  Eigen::LevenbergMarquardt<LMFunctor, double> lm(functor);
  // Set maximum number of function evaluations if it was set to a sensible value
  if (maxFuncEval > 0)
    lm.parameters.maxfev = maxFuncEval;
  lm.minimize(parameters);
  // Calculate and update covariance matrix if desired
  if (calculateCovarianceMatrix) {
    auto hessian = lm.fjac.transpose() * lm.fjac;
    auto inverseHessian = hessian.inverse();
    auto variance = (1.0 / (functor.m - functor.n + 1.0)) * lm.fvec.squaredNorm();
    covarianceMatrix_ = variance * inverseHessian;
  }
}

LevenbergMarquardt::LMFunctor::LMFunctor(UpdateFunctionManagerBase& updateFunctionManager)
  : updateFunctionManager_(updateFunctionManager) {
}

// Compute 'm' errors, one for each data point, for the given parameter values in 'x'
int LevenbergMarquardt::LMFunctor::operator()(const Eigen::VectorXd& parameters, Eigen::VectorXd& fvec) const {
  updateFunctionManager_.updateErrors(parameters, fvec);
  return 0;
}

// Compute the Jacobian of the errors
int LevenbergMarquardt::LMFunctor::df(const Eigen::VectorXd& parameters, Eigen::MatrixXd& fjac) const {
  updateFunctionManager_.updateJacobian(parameters, fjac);
  return 0;
}

// Returns 'm', the number of values.
int LevenbergMarquardt::LMFunctor::values() const {
  return m;
}

// Returns 'n', the number of inputs.
int LevenbergMarquardt::LMFunctor::inputs() const {
  return n;
}

const Eigen::MatrixXd& LevenbergMarquardt::getCovarianceMatrix() {
  return covarianceMatrix_;
}

} // namespace Utils
} // namespace Scine