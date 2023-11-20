/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ComplexSpinAdaptedMatrix.h"

namespace Scine {
namespace Utils {

ComplexSpinAdaptedMatrix ComplexSpinAdaptedMatrix::createRestricted(Eigen::MatrixXd m) {
  ComplexSpinAdaptedMatrix matrix;
  matrix.resize(static_cast<int>(m.rows()));
  matrix.setRestrictedMatrix(std::move(m));
  return matrix;
}

ComplexSpinAdaptedMatrix ComplexSpinAdaptedMatrix::createUnrestricted(Eigen::MatrixXd alpha, Eigen::MatrixXd beta) {
  assert(alpha.size() == beta.size() && "Alpha and beta matrix don't have the same size");
  ComplexSpinAdaptedMatrix matrix;
  matrix.resize(static_cast<int>(alpha.rows()));
  matrix.setAlphaMatrix(std::move(alpha));
  matrix.setBetaMatrix(std::move(beta));
  return matrix;
}

void ComplexSpinAdaptedMatrix::resize(int nAOs) {
  restrictedMatrix_.resize(nAOs, nAOs);
  alphaMatrix_.resize(nAOs, nAOs);
  betaMatrix_.resize(nAOs, nAOs);
}

} // namespace Utils
} // namespace Scine
