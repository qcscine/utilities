/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_LCAOUTIL_HFWAVEFUNCTIONOVERLAP_H
#define UTILS_LCAOUTIL_HFWAVEFUNCTIONOVERLAP_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {

class OccupiedMolecularOrbitals;

namespace LcaoUtil {

/*!
 * Class to calculate the overlap between two Hartree--Fock-like wave functions.
 */
class HFWaveFunctionOverlap {
 public:
  /*! Calculate the overlap in an orthonormal basis. */
  static double calculateOrthonormalOverlap(const OccupiedMolecularOrbitals& c1, const OccupiedMolecularOrbitals& c2);
  /*! Calculate the overlap in a non-orthonormal basis. */
  static double calculateNonOrthonormalOverlap(const OccupiedMolecularOrbitals& c1, const OccupiedMolecularOrbitals& c2,
                                               const Eigen::MatrixXd& s);

 private:
  static double unrestrictedOrthonormalOverlap(const OccupiedMolecularOrbitals& c1, const OccupiedMolecularOrbitals& c2);
  static double unrestrictedNonOrthonormalOverlap(const OccupiedMolecularOrbitals& c1,
                                                  const OccupiedMolecularOrbitals& c2, const Eigen::MatrixXd& s);
  static double orthonormalContribution(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2);
  static double nonOrthonormalContribution(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2, const Eigen::MatrixXd& s);
};

} // namespace LcaoUtil

} // namespace Utils
} // namespace Scine
#endif // LCAOUTIL_HFWAVEFUNCTIONOVERLAP_H