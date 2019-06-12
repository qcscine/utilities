/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_ENERGYWEIGHTEDDENSITYMATRIXBUILDER_H
#define UTILS_ENERGYWEIGHTEDDENSITYMATRIXBUILDER_H

#include <Eigen/Core>
#include <utility>
#include <vector>

namespace Scine {
namespace Utils {

class MolecularOrbitals;
class DensityMatrix;
class SingleParticleEnergies;
namespace LcaoUtil {

/*!
 * Class to generate energy-weighted density matrices for given occupations
 */
class EnergyWeightedDensityMatrixBuilder {
 public:
  explicit EnergyWeightedDensityMatrixBuilder(const MolecularOrbitals& coefficientMatrix,
                                              const SingleParticleEnergies& orbitalEnergies);

  /*! Generates a restricted energy-weighted density matrix for the given number of electrons. */
  Eigen::MatrixXd generateRestrictedForNumberElectrons(int nElectrons) const;
  /*! Generates an unrestricted energy-weighted density matrix for the given number of electrons and corresponding spin
   * multiplicity. */
  Eigen::MatrixXd generateUnrestrictedForNumberElectronsAndMultiplicity(int nElectrons, int spinMultiplicity) const;
  /*! Generates an unrestricted energy-weighted density matrix for the given number of alpha and beta electrons. */
  Eigen::MatrixXd generateUnrestrictedForNumberAlphaAndBetaElectrons(int nAlpha, int nBeta) const;
  /*! Generates a restricted energy-weighted density matrix from specified molecular orbitals (eigenvectors), which will
   * be doubly filled. */
  Eigen::MatrixXd generateRestrictedForSpecifiedOrbitals(const std::vector<int>& doublyOccupiedOrbitals) const;
  /*! Generates an unrestricted energy-weighted density matrix from specified molecular orbitals (eigenvectors). */
  Eigen::MatrixXd generateUnrestrictedForSpecifiedOrbitals(const std::vector<int>& alphaOrbitals,
                                                           const std::vector<int>& betaOrbitals) const;

 private:
  Eigen::MatrixXd calculateSingleOrbitalEWDensity(const Eigen::VectorXd& eigenvector, double orbitalEnergy) const;
  const MolecularOrbitals& coefficientMatrix_;
  const SingleParticleEnergies& orbitalEnergies_;
};

} // namespace LcaoUtil

} // namespace Utils
} // namespace Scine
#endif // UTILS_ENERGYWEIGHTEDDENSITYMATRIXBUILDER_H
