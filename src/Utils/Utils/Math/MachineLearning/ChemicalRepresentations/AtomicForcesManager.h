/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_ML_CHEMICALREPRESENTATIONS_ATOMICFORCESMANAGER_H
#define UTILS_ML_CHEMICALREPRESENTATIONS_ATOMICFORCESMANAGER_H

#include <Utils/Constants.h>
#include <Utils/Typenames.h>
#include <Eigen/Dense>

namespace Scine {
namespace Utils {
class AtomCollection;
namespace MachineLearning {

/**
 * @class AtomicForcesManager AtomicForcesManager.h
 * @brief This class manages the atomic forces of a molecular structure for machine learning applications.
 *        It constructs the transformation matrices to the internal force representations.
 *        Furthermore, it can calculate the features for a given atom for the machine learning of forces
 *        as described in arXiv:1908.10492.
 */
class AtomicForcesManager {
 public:
  /**
   * @brief Constructor. It constructs the transformation matrices to the molecular internal coordinate system.
   *        The transformation matrices are re-calculated each time a new set of positions is given.
   * @param structure The molecular structure.
   */
  explicit AtomicForcesManager(const AtomCollection& structure);
  /**
   * @brief Sets a new set of positions. It also updates the transformation matrices for all atoms.
   */
  void modifyPositions(const PositionCollection& newPositions);
  /**
   * @brief Transforms a force from an internal representation to a global representation for a given atom.
   */
  Eigen::RowVector3d toGlobalRepresentation(const Eigen::Ref<Eigen::RowVector3d> force, int atomIndex);
  /**
   * @brief Transforms a force from a global representation to an internal representation for a given atom.
   */
  Eigen::RowVector3d toInternalRepresentation(const Eigen::Ref<Eigen::RowVector3d> force, int atomIndex);
  /**
   * @brief Calculates and returns the feature vector for a given atom and the current
   *        set of positions of the molecular system.
   * @param atomIndex The index of the atom to return the feature vector for.
   * @return The elements of the feature matrix.
   */
  Eigen::VectorXd calculateFeatures(int atomIndex);

 private:
  // This function fills the importantAtoms_ vector based on the original structure
  // and is called in the constructor of this class.
  void determineImportantAtoms();
  // Calculates and returns the feature matrix for a given atom and the current set of positions of the molecular system.
  Eigen::MatrixXd calculateFeatureMatrix(int atomIndex);
  // Calculates and returns the center of nuclear charge of the molecular system in the global reference frame
  Eigen::RowVector3d calculateCenterOfCharge();
  // Constructs the transformation matrices for all atoms
  void constructTransformationMatrices();
  // Calculates and returns the chemical environment matrix for a given atom
  Eigen::MatrixXd calculateChemicalEnvironment(int atomIndex);
  // Sorts the rows of the feature matrix and returns the sorted matrix
  Eigen::MatrixXd sortFeatureMatrix(std::vector<std::pair<Eigen::RowVector3d, double>>& unsortedFeatures);
  // The transformation matrices for all atoms,
  // which convert the internal force representation to the global one
  std::vector<Eigen::MatrixXd> transformationMatrices_;
  // The molecular structure.
  const AtomCollection& structure_;
  // The current set of positions
  PositionCollection currentPositions_;
  // The number of atoms
  int nAtoms_;
  // A vector that holds a vector for each atom, which contains the atom indices of those atoms
  // that are closer than the distance threshold in the original structure.
  std::vector<std::vector<int>> importantAtoms_;
  // The distance threshold deciding which atoms are considered in the chemical environment matrix and the feature matrix
  static constexpr double distanceThreshold_ = 8.0 * Constants::bohr_per_angstrom;
};

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine

#endif // UTILS_ML_CHEMICALREPRESENTATIONS_ATOMICFORCESMANAGER_H
