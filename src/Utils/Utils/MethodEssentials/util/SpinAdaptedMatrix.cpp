/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SpinAdaptedMatrix.h"

namespace Scine {
namespace Utils {

SpinAdaptedMatrix SpinAdaptedMatrix::createRestricted(Eigen::MatrixXd m) {
  SpinAdaptedMatrix matrix;
  matrix.resize(static_cast<int>(m.rows()));
  matrix.setRestrictedMatrix(std::move(m));
  return matrix;
}

SpinAdaptedMatrix SpinAdaptedMatrix::createUnrestricted(Eigen::MatrixXd alpha, Eigen::MatrixXd beta) {
  assert(alpha.size() == beta.size() && "Alpha and beta matrix don't have the same size");
  SpinAdaptedMatrix matrix;
  matrix.resize(static_cast<int>(alpha.rows()));
  matrix.setAlphaMatrix(std::move(alpha));
  matrix.setBetaMatrix(std::move(beta));
  return matrix;
}

void SpinAdaptedMatrix::resize(int nAOs) {
  restrictedMatrix_.resize(nAOs, nAOs);
  alphaMatrix_.resize(nAOs, nAOs);
  betaMatrix_.resize(nAOs, nAOs);
}

} // namespace Utils
} // namespace Scine
