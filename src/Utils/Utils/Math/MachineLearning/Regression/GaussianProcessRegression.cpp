/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GaussianProcessRegression.h"
#include "RegressionSettings.h"
#include <Core/Log.h>
#include <LBFGSB.h>
#include <Utils/Constants.h>
#include <Utils/Optimizer/GradientBased/Bfgs.h>
#include <Utils/Optimizer/GradientBased/Lbfgs.h>
#include <Utils/Optimizer/GradientBased/SteepestDescent.h>
#include <Utils/Optimizer/LeastSquares/LevenbergMarquardt.h>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <random>

namespace Scine {
namespace Utils {
namespace MachineLearning {

GaussianProcessRegression::GaussianProcessRegression() {
  this->settings_ = std::make_unique<RegressionSettings>();
}

GaussianProcessRegression::GaussianProcessRegression(const GaussianProcessRegression& rhs) {
  auto valueCollection = dynamic_cast<const Utils::UniversalSettings::ValueCollection&>(rhs.settings());
  this->settings_ =
      std::make_unique<Utils::Settings>(Utils::Settings(valueCollection, rhs.settings().getDescriptorCollection()));
}

void GaussianProcessRegression::trainModel(const Eigen::MatrixXd& featureValues, const Eigen::MatrixXd& targetValues) {
  if (targetValues.rows() != featureValues.rows()) {
    throw std::runtime_error("The number of data points does not match between the feature and target matrices.");
  }

  if (targetValues.cols() > 1)
    throw std::runtime_error("Gaussian Process Regression is not implemented for n-dimensional target value input.");

  targetValues_ = targetValues.transpose();
  featureValues_ = featureValues.transpose();
  numDataPoints_ = static_cast<int>(targetValues_.cols());

  // Initialize Hyperparameters default
  hyperparameters_.parameters = Eigen::VectorXd::Constant(3, thetaSpecifier_.guess);
  hyperparameters_.parameters(1) = sigmaFSqSpecifier_.guess;
  hyperparameters_.parameters(2) = sigmaYSqSpecifier_.guess;
  // TODO: Set Hyperparameters from settings
  // Set Params to optimize
  hyperparameters_.toOptimize = Eigen::Matrix<bool, -1, 1>::Constant(3, true);
  hyperparameters_.toOptimize(0) = thetaSpecifier_.toOptimize;
  hyperparameters_.toOptimize(1) = sigmaFSqSpecifier_.toOptimize;
  hyperparameters_.toOptimize(2) = sigmaYSqSpecifier_.toOptimize;

  // maximise the log-marginal likelihood
  fit();
}

Eigen::VectorXd GaussianProcessRegression::predict(const Eigen::VectorXd& data) {
  numTestPoints_ = static_cast<int>(data.size());
  variance_ = Eigen::MatrixXd::Zero(numTestPoints_, numDataPoints_);
  // construct correlation matrices
  Eigen::MatrixXd kernelMatrix = constructCorrelationMatrix(featureValues_, featureValues_);
  kernelMatrix.diagonal().array() += hyperparameters_.parameters(2);
  // Ks
  Eigen::MatrixXd Ks = constructCorrelationMatrix(featureValues_, data.transpose());
  // Kss
  Eigen::MatrixXd Kss = constructCorrelationMatrix(data.transpose(), data.transpose());

  // compute mean of prediction
  // mu = Ks.T K_y^-1 * y
  Eigen::JacobiSVD<Eigen::MatrixXd> solver(kernelMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd y(Eigen::Map<const Eigen::VectorXd>(targetValues_.data(), targetValues_.size()));
  Eigen::VectorXd K_yInvY = solver.solve(y);
  Eigen::VectorXd meanOfPrediction = Ks.transpose() * K_yInvY;

  // predict the variance
  // Sigma = Kss - Ks.T * K_y^-1 Ks
  Eigen::MatrixXd secondTerm = Ks.transpose() * solver.solve(Ks); // K_yInvKs;
  variance_ = Eigen::MatrixXd::Zero(numTestPoints_, numTestPoints_);
  variance_ = Kss - secondTerm;
  return meanOfPrediction;
}

Eigen::MatrixXd GaussianProcessRegression::getVarianceOfPrediction() {
  return variance_;
}

void GaussianProcessRegression::fit() {
  // Optimize the hyperparameter of a gaussian kernel
  // external LBFGS-B optimizer
  LBFGSpp::LBFGSBParam<double> pm;
  pm.max_iterations = settings_->getInt(SettingsNames::Optimizations::MachineLearning::maxIterations);
  pm.max_linesearch = settings_->getInt(SettingsNames::Optimizations::MachineLearning::maxLinesearch);
  pm.min_step = 1e-20;
  pm.max_step = 1e+20;
  pm.epsilon = settings_->getDouble(SettingsNames::Optimizations::MachineLearning::convergenceTolerance);
  pm.ftol = settings_->getDouble(SettingsNames::Optimizations::MachineLearning::ftol);

  LBFGSpp::LBFGSBSolver<double> solver(pm);

  Eigen::VectorXd lb(3);
  Eigen::VectorXd ub(3);

  lb[0] = (thetaSpecifier_.toOptimize) ? thetaSpecifier_.bounds->first : thetaSpecifier_.guess;
  lb[1] = (sigmaFSqSpecifier_.toOptimize) ? sigmaFSqSpecifier_.bounds->first : sigmaFSqSpecifier_.guess;
  lb[2] = (sigmaYSqSpecifier_.toOptimize) ? sigmaYSqSpecifier_.bounds->first : sigmaYSqSpecifier_.guess;

  ub[0] = (thetaSpecifier_.toOptimize) ? thetaSpecifier_.bounds->second : thetaSpecifier_.guess;
  ub[1] = (sigmaFSqSpecifier_.toOptimize) ? sigmaFSqSpecifier_.bounds->second : sigmaFSqSpecifier_.guess;
  ub[2] = (sigmaYSqSpecifier_.toOptimize) ? sigmaYSqSpecifier_.bounds->second : sigmaYSqSpecifier_.guess;

  for (int i = 0; i < 3; i++) {
    if (lb[i] > ub[i])
      throw std::runtime_error("Upper bound cannot be smaller than lower bound!");
  }

  const auto logLikelihoodFunction = [&, this](const Eigen::VectorXd& parameters, Eigen::VectorXd& gradients) {
    // set loglikelihood value to 0
    double value = 0.0;
    // n parameters, n gradients
    const int nParams = hyperparameters_.nParams();
    gradients = Eigen::VectorXd::Zero(nParams);
    // set parameters value to the one given by the optimizer
    assert(parameters.size() == nParams);
    // With hard-coded isotropic gaussian kernel only 3 parameters.
    assert(nParams == 3);
    hyperparameters_.parameters = parameters;

    // Construct the kernel matrix and its derivative
    Eigen::MatrixXd kernelMatrix(numDataPoints_, numDataPoints_);
    std::vector<Eigen::MatrixXd> derivativeKernelMatrix{size_t(nParams), Eigen::MatrixXd::Zero(numDataPoints_, numDataPoints_)};
    AutomaticDifferentiation::FirstND sigmaYSq;
    if (sigmaYSqSpecifier_.toOptimize) {
      sigmaYSq = {hyperparameters_.parameters(nParams - 1), Eigen::Vector3d(0.0, 0.0, 1.0)};
    }
    else {
      sigmaYSq = {hyperparameters_.parameters(nParams - 1), Eigen::Vector3d(0.0, 0.0, 0.0)};
    }
    for (int i = 0; i < numDataPoints_; ++i) {
      for (int j = 0; j < numDataPoints_; ++j) {
        auto value = kernel_(featureValues_.col(i), featureValues_.col(j), hyperparameters_);
        if (i == j)
          value += sigmaYSq;
        kernelMatrix(i, j) = value.value();
        // Derivatives
        for (int param = 0; param < nParams; ++param) {
          // Add sigma_y derivative on the diagonal (deriv is either 0 or 1, depending on
          // whether sigma_y is being optimized)
          derivativeKernelMatrix[param](i, j) = value.derivative(param);
          // derivativeKernelMatrix[param](i, i) = value.derivative(param);
        }
      }
    }
    // kernelMatrix.diagonal().array() += hyperparameters_.parameters(nParams - 1);
    Eigen::JacobiSVD<Eigen::MatrixXd> solver(kernelMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double logDeterminant = 0.5 * std::log(kernelMatrix.fullPivLu().determinant());

    // compute negative log likelihood
    //
    // K_y^-1 * y, a row(K) * col(y) matrix. If y is 1-dim, then a nSamples * 1 matrix
    // K_y = K + sigma_y^2 * I
    Eigen::VectorXd y(Eigen::Map<const Eigen::VectorXd>(targetValues_.data(), targetValues_.size()));
    Eigen::VectorXd k_yInvY = solver.solve(y);
    value += 0.5 * y.dot(k_yInvY);
    value += 0.5 * logDeterminant;
    value += 0.5 * numDataPoints_ * log(2 * Utils::Constants::pi);

    // calculate derivative
    const auto calculateDerivative = [&](const Eigen::MatrixXd& derivativeMatrix) -> double {
      double logLikelihoodDer = 0.0;
      logLikelihoodDer = -k_yInvY.transpose() * derivativeMatrix * k_yInvY;
      logLikelihoodDer += 0.5 * solver.solve(derivativeMatrix).trace();
      return logLikelihoodDer;
    };

    for (int i = 0; i < nParams; ++i) {
      gradients(i) = calculateDerivative(derivativeKernelMatrix[i]);
    }

    return value;
  };

  double value = 0.0;
  solver.minimize(logLikelihoodFunction, hyperparameters_.parameters, value, lb, ub);
}

Eigen::VectorXd GaussianProcessRegression::getOptimizedHyperparameters() {
  return hyperparameters_.parameters;
}

void GaussianProcessRegression::setThetaHyperparametersGuess(HyperparameterSpecifier specifier) {
  thetaSpecifier_ = std::move(specifier);
}

void GaussianProcessRegression::setSigmaFSqHyperparametersGuess(HyperparameterSpecifier specifier) {
  sigmaFSqSpecifier_ = std::move(specifier);
}

void GaussianProcessRegression::setSigmaYSqHyperparametersGuess(HyperparameterSpecifier specifier) {
  sigmaYSqSpecifier_ = std::move(specifier);
}

auto GaussianProcessRegression::getBoundaries() const -> std::vector<std::pair<int, std::pair<double, double>>> {
  std::vector<std::pair<int, std::pair<double, double>>> boundaries;
  const int nParams = hyperparameters_.nParams();
  if (thetaSpecifier_.bounds) {
    for (int i = 0; i < nParams - 2; ++i) {
      boundaries.emplace_back(i, *thetaSpecifier_.bounds);
    }
  }
  if (sigmaFSqSpecifier_.bounds) {
    boundaries.emplace_back(nParams - 2, *sigmaFSqSpecifier_.bounds);
  }
  if (sigmaFSqSpecifier_.bounds) {
    boundaries.emplace_back(nParams - 1, *sigmaYSqSpecifier_.bounds);
  }

  return boundaries;
}

Eigen::MatrixXd GaussianProcessRegression::constructCorrelationMatrix(const Eigen::MatrixXd& data1,
                                                                      const Eigen::MatrixXd& data2) const {
  Eigen::MatrixXd K(data1.cols(), data2.cols());
  for (int i = 0; i < data1.cols(); ++i) {
    for (int j = 0; j < data2.cols(); ++j) {
      double value = kernel_(data1.col(i), data2.col(j), hyperparameters_).value();
      K(i, j) = value;
    }
  }
  return K;
}

const Utils::Settings& GaussianProcessRegression::settings() const {
  return *settings_;
}

Utils::Settings& GaussianProcessRegression::settings() {
  return *settings_;
}

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine
