/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_LCAOUTIL_H
#define UTILS_LCAOUTIL_H

#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {
class BondOrderCollection;
}
namespace Utils {

class SingleParticleEnergies;
class AtomsOrbitalsIndexes;
class DensityMatrix;
class SpinAdaptedMatrix;
class MolecularOrbitals;

/*! \file LcaoUtil.h
 * This header file contains the declarations for functions that are commonly used
 * by qc methods based on an LCAO ansatz.
 */

namespace LcaoUtil {

/*! Calculate the numbers of alpha and beta electrons from the total number of electrons and the spin multiplicity. */
void getNumberUnrestrictedElectrons(int& nAlpha, int& nBeta, int nElectrons, int spinMultiplicity);

void solveRestrictedEigenvalueProblem(const SpinAdaptedMatrix& fockMatrix, MolecularOrbitals& coefficientMatrix,
                                      SingleParticleEnergies& singleParticleEnergies);

void solveRestrictedGeneralizedEigenvalueProblem(const SpinAdaptedMatrix& fockMatrix, const Eigen::MatrixXd& overlapMatrix,
                                                 MolecularOrbitals& coefficientMatrix,
                                                 SingleParticleEnergies& singleParticleEnergies);

void solveUnrestrictedEigenvalueProblem(const SpinAdaptedMatrix& fockMatrix, MolecularOrbitals& coefficientMatrix,
                                        SingleParticleEnergies& singleParticleEnergies);

void solveUnrestrictedGeneralizedEigenvalueProblem(const SpinAdaptedMatrix& fockMatrix, const Eigen::MatrixXd& overlapMatrix,
                                                   MolecularOrbitals& coefficientMatrix,
                                                   SingleParticleEnergies& singleParticleEnergies);

void calculateRestrictedDensityMatrix(DensityMatrix& densityMatrix, const MolecularOrbitals& coefficientMatrix, int nElectrons);

void calculateUnrestrictedDensityMatrices(DensityMatrix& densityMatrix, const MolecularOrbitals& coefficientMatrix,
                                          int nElectrons, int spinMultiplicity);

void calculateRestrictedEnergyWeightedDensityMatrix(Eigen::MatrixXd& energyWeightedDensityMatrix,
                                                    const MolecularOrbitals& coefficientMatrix,
                                                    const SingleParticleEnergies& singleParticleEnergies, int nElectrons);

void calculateUnrestrictedEnergyWeightedDensityMatrix(Eigen::MatrixXd& energyWeightedDensityMatrix,
                                                      const MolecularOrbitals& coefficientMatrix,
                                                      const SingleParticleEnergies& singleParticleEnergies,
                                                      int nElectrons, int spinMultiplicity);

/*! Computes the lower-triangular bond-order matrix for an non-orthogonal basis. */
void calculateBondOrderMatrix(Utils::BondOrderCollection& bondOrderMatrix, const DensityMatrix& densityMatrix,
                              const Eigen::MatrixXd& overlapMatrix, const AtomsOrbitalsIndexes& aoIndexes);

/*! Computes the lower-triangular bond-order matrix for an orthogonal basis. */
void calculateOrthonormalBondOrderMatrix(Utils::BondOrderCollection& bondOrderMatrix,
                                         const DensityMatrix& densityMatrix, const AtomsOrbitalsIndexes& aoIndexes);

void calculateOrthonormalAtomicCharges(std::vector<double>& mullikenCharges, const std::vector<double>& coreCharges,
                                       const DensityMatrix& densityMatrix, const AtomsOrbitalsIndexes& aoIndexes);

void calculateMullikenAtomicCharges(std::vector<double>& mullikenCharges, const std::vector<double>& coreCharges,
                                    const DensityMatrix& densityMatrix, const Eigen::MatrixXd& overlapMatrix,
                                    const AtomsOrbitalsIndexes& aoIndexes);

} // namespace LcaoUtil

} // namespace Utils
} // namespace Scine
#endif // UTILS_LCAOUTIL_H