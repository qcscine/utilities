/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/OccupiedMolecularOrbitals.h>
#include <Utils/Scf/LcaoUtils/HFWaveFunctionOverlap.h>
#include <Eigen/Dense>

namespace Scine {
namespace Utils {

namespace LcaoUtils {

double HFWaveFunctionOverlap::calculateOrthonormalOverlap(const OccupiedMolecularOrbitals& c1,
                                                          const OccupiedMolecularOrbitals& c2) {
  // Convert to uhf to make comparison of RHF - UHF possible
  // TODO: make more efficient if no conversion is needed
  auto unrestricted1 = c1.toUnrestricted();
  auto unrestricted2 = c2.toUnrestricted();
  return unrestrictedOrthonormalOverlap(unrestricted1, unrestricted2);
}

double HFWaveFunctionOverlap::calculateNonOrthonormalOverlap(const OccupiedMolecularOrbitals& c1,
                                                             const OccupiedMolecularOrbitals& c2, const Eigen::MatrixXd& s) {
  // Convert to uhf to make comparison of RHF - UHF possible
  // TODO: make more efficient if no conversion is needed
  auto unrestricted1 = c1.toUnrestricted();
  auto unrestricted2 = c2.toUnrestricted();
  return unrestrictedNonOrthonormalOverlap(unrestricted1, unrestricted2, s);
}

double HFWaveFunctionOverlap::unrestrictedOrthonormalOverlap(const OccupiedMolecularOrbitals& c1,
                                                             const OccupiedMolecularOrbitals& c2) {
  double f1 = orthonormalContribution(c1.alphaMatrix(), c2.alphaMatrix());
  double f2 = orthonormalContribution(c1.betaMatrix(), c2.betaMatrix());
  // Return a product of determinants:
  // The determinant of a block diagonal matrix is the product of the determinants of the block
  return f1 * f2;
}

double HFWaveFunctionOverlap::unrestrictedNonOrthonormalOverlap(const OccupiedMolecularOrbitals& c1,
                                                                const OccupiedMolecularOrbitals& c2,
                                                                const Eigen::MatrixXd& s) {
  double f1 = nonOrthonormalContribution(c1.alphaMatrix(), c2.alphaMatrix(), s);
  double f2 = nonOrthonormalContribution(c1.betaMatrix(), c2.betaMatrix(), s);
  // Return a product of determinants:
  // The determinant of a block diagonal matrix is the product of the determinants of the block
  return f1 * f2;
}

double HFWaveFunctionOverlap::orthonormalContribution(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2) {
  if (m1.cols() != m2.cols()) {
    throw std::runtime_error("Not possible to calculate overlap between systems with different number of electrons.");
  }
  return (m1.transpose() * m2).determinant();
}

double HFWaveFunctionOverlap::nonOrthonormalContribution(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2,
                                                         const Eigen::MatrixXd& s) {
  if (m1.cols() != m2.cols()) {
    throw std::runtime_error("Not possible to calculate overlap between systems with different number of electrons.");
  }
  return (m1.transpose() * s * m2).determinant();
}

} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine
