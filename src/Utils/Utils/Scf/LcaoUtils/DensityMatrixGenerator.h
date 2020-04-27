/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_DENSITYMATRIXGENERATOR_H
#define UTILS_DENSITYMATRIXGENERATOR_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {

class MolecularOrbitals;
class DensityMatrix;
class SingleParticleEnergies;
namespace LcaoUtils {
class ElectronicOccupation;

/*!
 * Class to convert an electronic occupation and the corresponding molecular orbitals
 * to a density matrix.
 * \sa DensityMatrixBuilder
 */
class DensityMatrixGenerator {
 public:
  static DensityMatrix generate(const ElectronicOccupation& occupation, const MolecularOrbitals& molecularOrbitals);
  static Eigen::MatrixXd generateEnergyWeighted(const ElectronicOccupation& occupation,
                                                const MolecularOrbitals& molecularOrbitals,
                                                const SingleParticleEnergies& orbitalEnergies);

 private:
};

} // namespace LcaoUtils

} // namespace Utils
} // namespace Scine
#endif // UTILS_DENSITYMATRIXGENERATOR_H