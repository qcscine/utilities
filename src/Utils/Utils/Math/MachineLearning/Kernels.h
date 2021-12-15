/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_KERNELS_H
#define UTILS_MATH_KERNELS_H

#include <Eigen/Dense>
#include <vector>

namespace Scine {
namespace Utils {
namespace MachineLearning {

namespace Kernels {

/**
 * @brief The linear kernel. No hyperparameters.
 */
const std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&, const std::vector<double>&)> linearKernel =
    [](const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2, const std::vector<double>& /* hyperparameters */) {
      return vec1.dot(vec2);
    };

/**
 * @brief The Gaussian kernel. Only one hyperparameter.
 */
const std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&, const std::vector<double>&)> gaussianKernel =
    [](const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2, const std::vector<double>& hyperparameters) {
      return std::exp((-1 / (2 * std::pow(hyperparameters.at(0), 2))) * (vec1 - vec2).squaredNorm());
    };

/**
 * @brief The Laplacian kernel. Only one hyperparameter.
 */
const std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&, const std::vector<double>&)> laplacianKernel =
    [](const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2, const std::vector<double>& hyperparameters) {
      return std::exp((-1 / hyperparameters.at(0)) * (vec1 - vec2).lpNorm<1>());
    };

/**
 * @brief The polynomial kernel. Hyperparameters are, 0: degree of polynomial, 1: scaling factor for dot product,
 *                               2: constant added to scaled dot product.
 */
const std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&, const std::vector<double>&)> polynomialKernel =
    [](const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2, const std::vector<double>& hyperparameters) {
      return std::pow(hyperparameters.at(1) * vec1.dot(vec2) + hyperparameters.at(2), hyperparameters.at(0));
    };

} // namespace Kernels

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_KERNELS_H
