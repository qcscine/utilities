/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_SPINCONTAMINATION_H
#define UTILSOS_SPINCONTAMINATION_H

#include <Eigen/Dense>
#include <utility>

namespace Scine {
namespace Utils {
namespace LcaoUtils {

/*!
 * @brief Evaluates the spin contamination.
 *
 * @param overlap           Overlap matrix in AO basis.
 * @param coeffAlpha        MO coefficient matrix for alpha orbitals
 * @param coeffBeta         MO coefficient matrix for beta orbitals
 * @param nAlpha            Number of alpha electrons
 * @param nBeta             Number of beta electrons
 *
 * @return spin contamination and exact S^2 value
 */
inline static auto S2(const Eigen::MatrixXd& overlap, const Eigen::MatrixXd& coeffAlpha,
                      const Eigen::MatrixXd& coeffBeta, int nAlpha, int nBeta) -> std::pair<double, double> {
  auto Sab =
      coeffAlpha.block(0, 0, overlap.rows(), nAlpha).transpose() * overlap * coeffBeta.block(0, 0, overlap.rows(), nBeta);

  double exact = (nAlpha - nBeta) / 2. * ((nAlpha - nBeta) / 2. + 1);

  return {exact + nBeta - Sab.array().square().sum(), exact};
}

} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine

#endif // UTILSOS_SPINCONTAMINATION_H
