/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_MOLECULARTRAJECTORY_H
#define UTILS_MOLECULARTRAJECTORY_H

#include "Utils/DataStructures/PeriodicBoundaries.h"
#include "Utils/Typenames.h"

namespace Scine {
namespace Utils {

/**
 * @class MolecularTrajectory MolecularTrajectory.h
 * @brief A container class for Scine::Utils::PositionCollection classes.
 *
 * Mainly an override of std::vector<PositionCollection>, plus Element Type information.
 * Generally used for storing the positions of several consecutive molecular structures.
 */
class MolecularTrajectory {
 public:
  using Container = std::vector<PositionCollection>;
  // The energy container is always either kept empty or the same size as the container of structures.
  using EnergyContainer = std::vector<double>;
  using PbcContainer = std::vector<Eigen::Matrix3d>;
  using iterator = Container::iterator;
  using const_iterator = Container::const_iterator;
  using reference = Container::reference;
  using const_reference = Container::const_reference;

  /// @brief Default constructor.
  MolecularTrajectory() = default;
  explicit MolecularTrajectory(const ElementTypeCollection& elements);
  explicit MolecularTrajectory(double minimumRmsdForAddition);
  explicit MolecularTrajectory(const ElementTypeCollection& elements, double minimumRmsdForAddition);
  /**
   * @brief Set a single element type.
   * @param index The index of the atom to be altered
   * @param e The element type to be set.
   */
  void setElementType(int index, ElementType e);
  /**
   * @brief Set the all element types.
   * @param elements The new element types.
   */
  void setElementTypes(const ElementTypeCollection& elements);
  /**
   * @brief Get all the element types.
   * @return const ElementTypeCollection& A reference to all element types.
   */
  const ElementTypeCollection& getElementTypes() const;
  /**
   * @brief Setter for the energies.
   *
   *        Energies can only be set if the size of the energy vector
   *        matches the size of the structure vector.
   *
   * @param energies The energies.
   */
  void setEnergies(const EnergyContainer& energies);
  /**
   * @brief Setter for the periodic boundaries.
   *
   *        Pbcs can only be set if the size of the pbc vector
   *        matches the size of the structure vector.
   *
   * @param pbcs The pbcs in a vector as matrices.
   */
  void setPbcs(const PbcContainer& pbcs);
  /**
   * @brief Setter for the periodic boundaries.
   *
   *        Pbcs can only be set if the size of the pbc vector
   *        matches the size of the structure vector.
   *
   * @param pbcs The pbcs in a vector as true PeriodicBoundaries objects.
   */
  void setPbcs(const std::vector<PeriodicBoundaries>& pbcs);
  /**
   * @brief Getter for the vector of energies corresponding to the structures.
   * @return The vector of energies.
   */
  EnergyContainer getEnergies() const;
  /**
   * @brief Getter for the vector of periodic boundaries corresponding to the structures.
   * @return The vector of pbcs.
   */
  std::vector<PeriodicBoundaries> getPbcs() const;
  /**
   * @brief Clears the energy container.
   */
  void clearEnergies();
  /**
   * @brief Clears the pbc container.
   */
  void clearPbcs();
  /**
   * @brief Clears all steps in the trajectory, element information is retained.
   */
  void clear();
  /**
   * @brief Resets the number of PositionCollections.
   * @param n Number of PositionCollections.
   */
  void resize(int n);
  /**
   * @brief Adds a new set of positions (a frame)
   *
   *        This function is only usable if no energies have been previously added,
   *        thus the energy vector has to be empty.
   *
   * @param p The new positions to be added.
   */
  void push_back(PositionCollection p);
  /**
   * @brief Adds a new set of positions (a frame) with a corresponding energy.
   *
   *        This function is only usable if the energy vector already
   *        has the same size as the structure vector. Therefore, if a few
   *        structures are added without an energy, the energies have to be
   *        first set for the other structures, before this function becomes
   *        applicable.
   *
   * @param p The new positions to be added.
   * @param e The corresponding energy.
   */
  void push_back(PositionCollection p, double e);
  /**
   * @brief Adds a new set of positions (a frame) with a corresponding pbc.
   *
   *        This function is only usable if the pbc vector already
   *        has the same size as the structure vector. Therefore, if a few
   *        structures are added without a pbc, the pbcs have to be
   *        first set for the other structures, before this function becomes
   *        applicable.
   *
   * @param p The new positions to be added.
   * @param pbc The corresponding pbc.
   */
  void push_back(PositionCollection p, PeriodicBoundaries pbc);
  /**
   * @brief Adds a new set of positions (a frame) with a corresponding energy
   *        and periodic boundaries.
   *
   *        This function is only usable if the energy and pbc vector already
   *        have the same size as the structure vector. Therefore, if a few
   *        structures are added without an energy or pbc, they have to be
   *        first set for the other structures, before this function becomes
   *        applicable.
   *
   * @param p The new positions to be added.
   * @param e The corresponding energy.
   * @param pbc The corresponding periodic boundaries.
   */
  void push_back(PositionCollection p, double e, PeriodicBoundaries pbc);
  /**
   * @brief Checks if there are structures present.
   * @return true If no PositionCollections are stored.
   * @return false If any PositionCollections are present.
   */
  bool empty() const;
  /**
   * @brief Getter for the number of structures.
   * @return int Returns the number of stored structures/frames/PositionCollections.
   */
  int size() const;
  /**
   * @brief Getter for the size of the structure.
   * @return int Returns number of atoms in the structure.
   */
  int molecularSize() const;
  /**
   * @brief An iterator pointing to the start of the structures (not the elements!).
   * @return iterator The iterator.
   */
  iterator begin();
  /**
   * @brief The const version of the start() iterator.
   * @return const_iterator The iterator.
   */
  const_iterator begin() const;
  /**
   * @brief Erases one PositionCollection from the stored list.
   * @param position The position to erase.
   * @return iterator An iterator to the nex position in the new list.
   */
  iterator erase(iterator position);
  /**
   * @brief An iterator pointing to the end of the structures (not the elements!).
   * @return iterator The iterator.
   */
  iterator end();
  /**
   * @brief The const version of the end() iterator.
   * @return const_iterator The iterator.
   */
  const_iterator end() const;
  /**
   * @brief The access operator.
   * @param i The index to be accessed.
   * @return reference A reference to the structure at index i.
   */
  reference operator[](int i);
  /**
   * @brief The const version of the access operator.
   * @param i The index to be accessed.
   * @return const_reference A const reference to the structure at index i.
   */
  const_reference operator[](int i) const;
  /**
   * @brief The at() function to access single structures in the trajectory.
   * @param i The index to be accessed.
   * @return reference A reference to the structure at index i.
   */
  reference at(int i);
  /**
   * @brief The const version of the at() function.
   * @param i The index to be accessed.
   * @return const_reference A const reference to the structure at index i.
   */
  const_reference at(int i) const;
  /**
   * @brief A function referencing the front of the list of structure.
   * @return reference A reference to the first entry in the list of structures.
   */
  reference front();
  /**
   * @brief The const version of the front() function.
   * @return const_reference A const reference to the first entry in the list of structures.
   */
  const_reference front() const;
  /**
   * @brief A function referencing the back of the list of structure.
   * @return reference A reference to the last entry in the list of structures.
   */
  reference back();
  /**
   * @brief The const version of the back() function.
   * @return const_reference A const reference to the last entry in the list of structures.
   */
  const_reference back() const;
  /**
   * @brief Multiplication assignment operator (for example for unit conversion)
   * @param f The scalar to multiply with.
   * @return const MolecularTrajectory& The trajectory, in-place.
   */
  const MolecularTrajectory& operator*=(double f);
  /**
   * @brief Division assignment operator (for example for unit conversion)
   * @param f The scalar to divide by.
   * @return const MolecularTrajectory& The trajectory, in-place.
   */
  const MolecularTrajectory& operator/=(double f);
  /**
   * @brief Multiplication operator (for example for unit conversion)
   * @param f The scalar to multiply with.
   * @return MolecularTrajectory The resulting trajectory.
   */
  MolecularTrajectory operator*(double f) const;
  /**
   * @brief Division operator (for example for unit conversion)
   * @param f The scalar to divide by.
   * @return MolecularTrajectory The resulting trajectory.
   */
  MolecularTrajectory operator/(double f) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html

      private :
    /**
     * @brief Checks if a reset of the elements with a given collection is allowed.
     * @param ec The Element collection.
     * @return true If the reset is allowed.
     * @return false If the reset is not allowed.
     */
    bool
    resettingElementTypeCollectionIsAllowed(const ElementTypeCollection& ec) const;
  /**
   * @brief Checks if the addition of a given collection of positions is allowed.
   * @param p The collection of positions.
   * @return true If the addition is allowed.
   * @return false If the addition is not allowed.
   */
  bool additionOfPositionCollectionIsAllowed(const PositionCollection& p) const;
  bool additionOfPositionCollectionIsAllowedBasedOnRmsd(const PositionCollection& p) const;
  Container structureVector_;
  ElementTypeCollection elements_;
  EnergyContainer energies_;
  PbcContainer pbcs_;
  double minMeanSquareDeviation_;
  bool respectMinRmsd_ = false;
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_MOLECULARTRAJECTORY_H
