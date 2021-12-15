/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_ML_CHEMICALREPRESENTATIONS_COULOMBMATRIX_H
#define UTILS_ML_CHEMICALREPRESENTATIONS_COULOMBMATRIX_H

#include <Eigen/Dense>

namespace Scine {
namespace Utils {
class AtomCollection;
namespace MachineLearning {

/**
 * @class CoulombMatrix CoulombMatrix.h
 * @brief Constructs and stores a Coulomb matrix for application in machine learning.
 *
 *        Reference: M. Rupp et al., Phys. Rev. Lett., 2012, 108, 058301.
 *
 *        Note that this class neither performs padding of the feature vector with zeros
 *        nor sorting of the Coulomb matrix in any way (e.g., by the bag of bonds strategy)
 *        and therefore is tailored towards usage in machine-learned force fields, where only
 *        one type of molecule is considered. To treat more than one molecule in an ML project,
 *        one needs to extend the functionalities of this class.
 *
 */
class CoulombMatrix {
 public:
  /**
   * @brief Constructor. It constructs a Coulomb matrix for a given molecular structure.
   * @param structure The molecular structure.
   */
  explicit CoulombMatrix(const AtomCollection& structure);
  /**
   * @brief Getter for the Coulomb matrix.
   */
  Eigen::MatrixXd getMatrix() const;
  /**
   * @brief Getter for the feature vector corresponding to the elements of the upper triangle of the Coulomb matrix.
   */
  Eigen::VectorXd getFeatures() const;

 private:
  // Constructs the Coulomb matrix (called in the constructor of this class)
  static Eigen::VectorXd generateCoulombMatrix(const AtomCollection& structure);
  // Returns the distance between the two atoms with indices i and j
  static double interatomicDistance(int i, int j, const AtomCollection& structure);
  // The feature vector corresponding to the elements of the upper triangle of the Coulomb matrix
  Eigen::VectorXd featureVector_;
  // The number of atoms
  int nAtoms_;
  // The exponent in the term for the diagonal elements of the Coulomb matrix
  static constexpr double diagonalElementsExponent_ = 2.4;
};

} // namespace MachineLearning
} // namespace Utils
} // namespace Scine

#endif // UTILS_ML_CHEMICALREPRESENTATIONS_COULOMBMATRIX_H
