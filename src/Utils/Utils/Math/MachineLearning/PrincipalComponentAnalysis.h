/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_PRINCIPALCOMPONENTANALYSIS_H
#define UTILS_MATH_PRINCIPALCOMPONENTANALYSIS_H

#include <Eigen/Dense>

namespace Scine {
namespace Utils {
namespace MachineLearning {

/// @brief Alias for a pair of principal components and their respective proportions of the total explained variance.
using PcaContainer = std::pair<Eigen::MatrixXd, Eigen::VectorXd>;

/**
 * @class PrincipalComponentAnalysis PrincipalComponentAnalysis.h
 * @brief This class performs a principal component analysis (PCA) applying the self-adjoint eigensolver from Eigen.
 */
class PrincipalComponentAnalysis {
 public:
  /**
   * @brief Constructor.
   * @param data The data on which to perfom the PCA. The number of rows is equal to the number of samples
   *             and the number of columns to the number of features.
   */
  explicit PrincipalComponentAnalysis(const Eigen::MatrixXd& data);

  /**
   * @brief Performs the PCA.
   * @param numComponents Number of components to consider in the PCA.
   * @return A std::pair containing the first numComponents principal components and their
   *         respective proportions of the total explained variance.
   *
   * Note that the returned components and proportions are sorted in an descending order from left to right by
   * the size of the proportions of explained variance.
   *
   */
  PcaContainer calculate(int numComponents) const;

 private:
  // The data matrix.
  const Eigen::MatrixXd& data_;
};

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_PRINCIPALCOMPONENTANALYSIS_H
