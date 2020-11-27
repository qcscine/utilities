/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "PrincipalComponentAnalysis.h"

namespace Scine {
namespace Utils {
namespace MachineLearning {

PrincipalComponentAnalysis::PrincipalComponentAnalysis(const Eigen::MatrixXd& data) : data_(data) {
  assert(data_.size() != 0);
}

PcaContainer PrincipalComponentAnalysis::calculate(int numComponents) const {
  Eigen::MatrixXd centeredData = data_.rowwise() - data_.colwise().mean();
  Eigen::MatrixXd covarianceMatrix = centeredData.adjoint() * centeredData;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(covarianceMatrix);

  // Return components and the proportions of the total variance which they explain
  Eigen::MatrixXd components = eig.eigenvectors().rightCols(numComponents);
  Eigen::VectorXd explainedVarianceRatios = eig.eigenvalues().tail(numComponents) / eig.eigenvalues().sum();
  return {components.rowwise().reverse(), (explainedVarianceRatios.reverse())};
}

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine