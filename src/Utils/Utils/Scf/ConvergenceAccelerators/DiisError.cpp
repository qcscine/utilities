/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DiisError.h"
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>

namespace Scine {
namespace Utils {

void DiisError::resize(unsigned subspaceSize) {
  errorMatrices.resize(subspaceSize);
}

void DiisError::setErrorFromMatrices(int index, const SpinAdaptedMatrix& fock, const DensityMatrix& density,
                                     const Eigen::MatrixXd& overlap) {
  if (unrestricted_) {
    errorMatrices[index] = calculateUnrestrictedErrorMatrix(fock, density, overlap);
  }
  else {
    errorMatrices[index] = calculateRestrictedErrorMatrix(fock, density, overlap);
  }
}

Eigen::MatrixXd DiisError::calculateUnrestrictedErrorMatrix(const SpinAdaptedMatrix& fock, const DensityMatrix& density,
                                                            const Eigen::MatrixXd& overlap) const {
  if (orthogonal_) {
    Eigen::MatrixXd alphaE = calculateOrthogonalErrorMatrix(fock.alphaMatrix(), density.alphaMatrix());
    Eigen::MatrixXd betaE = calculateOrthogonalErrorMatrix(fock.betaMatrix(), density.betaMatrix());
    return alphaE + betaE;
  }

  Eigen::MatrixXd alphaE = calculateErrorMatrix(fock.alphaMatrix(), overlap, density.alphaMatrix());
  Eigen::MatrixXd betaE = calculateErrorMatrix(fock.betaMatrix(), overlap, density.betaMatrix());
  return alphaE + betaE;
}

Eigen::MatrixXd DiisError::calculateRestrictedErrorMatrix(const SpinAdaptedMatrix& fock, const DensityMatrix& density,
                                                          const Eigen::MatrixXd& overlap) const {
  if (orthogonal_) {
    Eigen::MatrixXd restrictedE = calculateOrthogonalErrorMatrix(fock.restrictedMatrix(), density.restrictedMatrix());
    return restrictedE;
  }

  Eigen::MatrixXd restrictedE = calculateErrorMatrix(fock.restrictedMatrix(), overlap, density.restrictedMatrix());
  return restrictedE;
}

Eigen::MatrixXd DiisError::calculateErrorMatrix(const Eigen::MatrixXd& fock, const Eigen::MatrixXd& overlap,
                                                const Eigen::MatrixXd& density) {
  return fock.selfadjointView<Eigen::Lower>() * density * overlap - overlap * density * fock.selfadjointView<Eigen::Lower>();
}

Eigen::MatrixXd DiisError::calculateOrthogonalErrorMatrix(const Eigen::MatrixXd& fock, const Eigen::MatrixXd& density) {
  return fock.selfadjointView<Eigen::Lower>() * density - density * fock.selfadjointView<Eigen::Lower>();
}

double DiisError::getError(int i1, int i2) const {
  const auto& e1 = errorMatrices[i1];
  const auto& e2 = errorMatrices[i2];
  return e1.cwiseProduct(e2).sum();
}
} // namespace Utils
} // namespace Scine
