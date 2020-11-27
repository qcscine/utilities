/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_DENSITYMATRIXBUILDER_H
#define UTILS_DENSITYMATRIXBUILDER_H

#include <Utils/Scf/LcaoUtils/MolecularOrbitalsManipulation.h>
#include <Eigen/Core>
#include <utility>
#include <vector>

namespace Scine {
namespace Utils {

class MolecularOrbitals;
class DensityMatrix;
namespace LcaoUtils {

/*!
 * Class to generate density matrices from coefficient matrices.
 * Different possibilities to do so:
 * * Specify the number of electrons
 * * Specify which orbitals to consider
 * * Specify the number of electrons and some number of swaps
 * TODO: Use the class MolecularOrbitalsManipulation instead of reimplementing the mixes and swaps
 */
class DensityMatrixBuilder {
 public:
  /*! Alias for partly occupied orbital, given by orbital index and its occupation [0-1] (UHF) vs [0-2] (RHF). */
  using PartlyOccupiedOrbital = std::pair<int, double>;

  explicit DensityMatrixBuilder(const MolecularOrbitals& coefficientMatrix);
  /*! Generates a restricted density matrix for the given number of electrons. */
  DensityMatrix generateRestrictedForNumberElectrons(int nElectrons) const;
  /*! Generates an unrestricted density matrix for the given number of electrons and corresponding spin multiplicity. */
  DensityMatrix generateUnrestrictedForNumberElectronsAndMultiplicity(int nElectrons, int spinMultiplicity) const;
  /*! Generates an unrestricted density matrix for the given number of alpha and beta electrons. */
  DensityMatrix generateUnrestrictedForNumberAlphaAndBetaElectrons(int nAlpha, int nBeta) const;
  /*! Generates a restricted density matrix from specified molecular orbitals (eigenvectors), which will be doubly
   * filled. */
  DensityMatrix generateRestrictedForSpecifiedOrbitals(const std::vector<int>& doublyOccupiedOrbitals) const;
  /*! Generates an unrestricted density matrix from specified molecular orbitals (eigenvectors). */
  DensityMatrix generateUnrestrictedForSpecifiedOrbitals(const std::vector<int>& alphaOrbitals,
                                                         const std::vector<int>& betaOrbitals) const;
  /*! Generates a restricted density matrix with given Swaps.
   * \pre nElectrons is even. */
  DensityMatrix generateRestrictedWithSwaps(const std::vector<MolecularOrbitalsManipulation::DeprecatedSwap>& swaps,
                                            int nElectrons) const;
  /*! Generates an unrestricted density matrix with given Swaps. */
  DensityMatrix generateUnrestrictedWithSwaps(const std::vector<MolecularOrbitalsManipulation::DeprecatedSwap>& alphaSwaps,
                                              const std::vector<MolecularOrbitalsManipulation::DeprecatedSwap>& betaSwaps,
                                              int nAlpha, int nBeta) const;
  /*! Generates a restricted density matrix from specified molecular orbitals (eigenvectors) and their occupation [0-2].
   */
  DensityMatrix generateRestrictedForSpecifiedPartlyOccupiedOrbitals(const std::vector<PartlyOccupiedOrbital>& orbitals);
  /*! Generates an unrestricted density matrix from specified molecular orbitals (eigenvectors) and their occupation
   * [0-1]. */
  DensityMatrix generateUnrestrictedForSpecifiedPartlyOccupiedOrbitals(const std::vector<PartlyOccupiedOrbital>& alphaOrbitals,
                                                                       const std::vector<PartlyOccupiedOrbital>& betaOrbitals);
  /*! Generates a restricted density matrix with given Mixing.
   * \pre nElectrons is even. */
  DensityMatrix generateRestrictedWithMixing(const std::vector<MolecularOrbitalsManipulation::DeprecatedMix>& mix,
                                             int nElectrons) const;
  /*! Generates an unrestricted density matrix with given Mixing. */
  DensityMatrix generateUnrestrictedWithMixing(const std::vector<MolecularOrbitalsManipulation::DeprecatedMix>& alphaMix,
                                               const std::vector<MolecularOrbitalsManipulation::DeprecatedMix>& betaMix,
                                               int nAlpha, int nBeta) const;

 private:
  Eigen::MatrixXd calculateSingleOrbitalDensity(const Eigen::VectorXd& eigenvector) const;
  Eigen::MatrixXd calculateBlockOrbitalDensity(const Eigen::MatrixXd& eigenvectors) const;
  Eigen::MatrixXd calculateDensityMatrix(const Eigen::MatrixXd& coefficientMatrix, int nOccupiedLevels) const;
  Eigen::MatrixXd calculateDifferenceSwapDensity(const Eigen::MatrixXd& coefficientMatrix,
                                                 const std::vector<MolecularOrbitalsManipulation::DeprecatedSwap>& swaps,
                                                 int homoIndex) const;
  Eigen::MatrixXd calculateDifferenceMixDensity(const Eigen::MatrixXd& coefficientMatrix,
                                                const std::vector<MolecularOrbitalsManipulation::DeprecatedMix>& mix,
                                                int homoIndex) const;
  const MolecularOrbitals& coefficientMatrix_;
};

} // namespace LcaoUtils

} // namespace Utils
} // namespace Scine
#endif // UTILS_DENSITYMATRIXBUILDER_H