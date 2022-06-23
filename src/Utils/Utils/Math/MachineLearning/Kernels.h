/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_KERNELS_H
#define UTILS_MATH_KERNELS_H

#include "Utils/Math/AutomaticDifferentiation/FirstND.h"
#include "Utils/Math/DerivOrderEnum.h"
#include <Eigen/Dense>
#include <vector>

namespace Scine {
namespace Utils {
namespace MachineLearning {

namespace Kernels {

const std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&)> linearKernel =
    [](const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2, const Eigen::VectorXd& /* hyperparameters */) {
      return vec1.dot(vec2);
    };

/**
 * @brief The Gaussian kernel. Only one hyperparameter.
 */
const std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&)> gaussianKernel =
    [](const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2, const Eigen::VectorXd& hyperparameters) {
      return std::exp((-1 / (2 * std::pow(hyperparameters(0), 2))) * (vec1 - vec2).squaredNorm());
    };

struct Hyperparameters {
  auto generateDerivable() const -> std::vector<AutomaticDifferentiation::FirstND> {
    assert(toOptimize.size() == parameters.size());

    std::vector<AutomaticDifferentiation::FirstND> result;
    for (int i = 0; i < nParams(); ++i) {
      if (toOptimize(i)) {
        result.emplace_back(parameters(i), Eigen::VectorXd(Eigen::MatrixXd::Identity(nParams(), nParams()).row(i)));
      }
      else {
        result.emplace_back(parameters(i), Eigen::VectorXd::Zero(nParams()));
      }
    }
    return result;
  }
  auto nParams() const -> int {
    return int(parameters.size());
  }
  Eigen::VectorXd parameters;
  Eigen::Matrix<bool, -1, 1> toOptimize;
};
/**
 * @brief The isotropic Gaussian kernel. Only 3 hyperparameters.
 */
const std::function<Utils::AutomaticDifferentiation::FirstND(const Eigen::VectorXd&, const Eigen::VectorXd&, const Hyperparameters&)> autoGaussianKernel =
    [](const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2, const Hyperparameters& hyperparameters) {
      using namespace Utils::AutomaticDifferentiation;
      // Generate derivable parameters
      auto derivableParams = hyperparameters.generateDerivable();
      assert(derivableParams.size() == hyperparameters.nParams());
      // Get kernel,
      // derivParams[0] = theta
      // derivParams[1] = sigma_f^2
      // derivParams[2] = sigma_y^2
      FirstND kernel = square(derivableParams[1]) * exp(-0.5 * (vec1 - vec2).squaredNorm() / square(derivableParams[0]));
      return kernel;
    };

/**
 * @brief The anysotropic Gaussian kernel. 1 hyperparameter per feature + 2 hyperparameters.
 */
const std::function<Utils::AutomaticDifferentiation::FirstND(const Eigen::VectorXd&, const Eigen::VectorXd&, const Hyperparameters&)> autoAnisotropicGaussianKernel =
    [](const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2, const Hyperparameters& hyperparameters) {
      using namespace Utils::AutomaticDifferentiation;
      // Generate derivable parameters
      auto derivableParams = hyperparameters.generateDerivable();
      assert(derivableParams.size() == hyperparameters.nParams());
      // Calculate distance with feature-dependent scale length
      FirstND distance{0, Eigen::VectorXd::Zero(hyperparameters.nParams())};
      for (int i = 0; i < vec1.rows(); ++i) {
        distance += std::pow(vec1(i) - vec2(i), 2) / derivableParams[i];
      }
      // Get kernel,
      // derivParams[hyperparameters.nParams() - 2] = sigma_f^2
      // derivParams[hyperparameters.nParams() - 1] = sigma_y^2
      FirstND kernel = derivableParams[hyperparameters.nParams() - 2] * exp(-0.5 * distance) +
                       derivableParams[hyperparameters.nParams() - 1];
      return kernel;
    };
/**
 * @brief The Laplacian kernel. Only one hyperparameter.
 */
const std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&)> laplacianKernel =
    [](const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2, const Eigen::VectorXd& hyperparameters) {
      return std::exp((-1 / hyperparameters(0)) * (vec1 - vec2).lpNorm<1>());
    };

/**
 * @brief The polynomial kernel. Hyperparameters are, 0: degree of polynomial, 1: scaling factor for dot product,
 *                               2: constant added to scaled dot product.
 */
const std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&)> polynomialKernel =
    [](const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2, const Eigen::VectorXd& hyperparameters) {
      return std::pow(hyperparameters(1) * vec1.dot(vec2) + hyperparameters(2), hyperparameters(0));
    };

} // namespace Kernels

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_KERNELS_H
