/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SpinAdaptedVector.h"

namespace Scine {
namespace Utils {

SpinAdaptedVector SpinAdaptedVector::createRestricted(Eigen::VectorXd m) {
  SpinAdaptedVector matrix;
  matrix.resize(static_cast<int>(m.rows()));
  matrix.setRestrictedVector(std::move(m));
  return matrix;
}

SpinAdaptedVector SpinAdaptedVector::createUnrestricted(Eigen::VectorXd alpha, Eigen::VectorXd beta) {
  assert(alpha.size() == beta.size() && "Alpha and beta matrix don't have the same size");
  SpinAdaptedVector matrix;
  matrix.resize(static_cast<int>(alpha.rows()));
  matrix.setAlphaVector(std::move(alpha));
  matrix.setBetaVector(std::move(beta));
  return matrix;
}

void SpinAdaptedVector::resize(int nAOs) {
  restrictedVector_.resize(nAOs);
  alphaVector_.resize(nAOs);
  betaVector_.resize(nAOs);
}

} // namespace Utils
} // namespace Scine
