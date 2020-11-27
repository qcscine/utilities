/**
 * @file
 * @brief A file containing definitions of classes that are just different names
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_OCCUPIEDMOLECULARORBITALS_H
#define UTILS_OCCUPIEDMOLECULARORBITALS_H

#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <vector>

namespace Scine {
namespace Utils {

namespace LcaoUtils {
class ElectronicOccupation;
}
class MolecularOrbitals;

/*!
 * Class for the occupied part of the coefficient (molecular orbital) matrix.
 * Contains only the occupied eigenfunctions (molecular orbitals).
 * The dimension of the matrices is the number of basis functions (rows) and the number of occupied orbitals (columns)
 * TODO: Take SingleParticleEnergies inside this class?
 * \sa MolecularOrbitals
 */
class OccupiedMolecularOrbitals {
 public:
  using Matrix = SpinAdaptedMatrix::Matrix;
  OccupiedMolecularOrbitals() = default;
  /*! Construct from all the molecular orbitals and the electronic occupation. */
  OccupiedMolecularOrbitals(const MolecularOrbitals& allOrbitals, const LcaoUtils::ElectronicOccupation& occupation);

  /*! If the molecular orbitals are restricted, transforms them into unrestricted ones. */
  void makeUnrestricted();
  /*! Return a copy of the orbitals, transformed to the unrestricted variant if needed. */
  OccupiedMolecularOrbitals toUnrestricted() const;

  bool isRestricted() const;
  bool isUnrestricted() const;

  const Matrix& restrictedMatrix() const;
  const Matrix& alphaMatrix() const;
  const Matrix& betaMatrix() const;

 private:
  void constructRestricted(const MolecularOrbitals& allOrbitals, const LcaoUtils::ElectronicOccupation& occupation);
  void constructUnrestricted(const MolecularOrbitals& allOrbitals, const LcaoUtils::ElectronicOccupation& occupation);
  static Eigen::MatrixXd calculateMatrixForFilledOrbitals(const Eigen::MatrixXd& matrixWithAllOrbitals,
                                                          const std::vector<int>& filledOrbitalsIndexes);

  SpinAdaptedMatrix matrix_;
  bool unrestricted_ = false;
};

inline bool OccupiedMolecularOrbitals::isUnrestricted() const {
  return unrestricted_;
}

inline bool OccupiedMolecularOrbitals::isRestricted() const {
  return !isUnrestricted();
}

inline const OccupiedMolecularOrbitals::Matrix& OccupiedMolecularOrbitals::restrictedMatrix() const {
  return matrix_.restrictedMatrix();
}

inline const OccupiedMolecularOrbitals::Matrix& OccupiedMolecularOrbitals::alphaMatrix() const {
  return matrix_.alphaMatrix();
}

inline const OccupiedMolecularOrbitals::Matrix& OccupiedMolecularOrbitals::betaMatrix() const {
  return matrix_.betaMatrix();
}

} // namespace Utils
} // namespace Scine
#endif // OCCUPIEDMOLECULARORBITALS_H