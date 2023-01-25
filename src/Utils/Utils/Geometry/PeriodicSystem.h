/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_PERIODICSYSTEM_H_
#define UTILS_PERIODICSYSTEM_H_

#include "Utils/DataStructures/PeriodicBoundaries.h"
#include "Utils/Geometry/Atom.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Typenames.h"
#include <Core/BaseClasses/ObjectWithLog.h>
#include <Eigen/Core>
#include <functional>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace Scine {
namespace Utils {

/**
 * @class PeriodicSystem PeriodicSystem.h
 * @brief A collection of Atoms including periodic boundary conditions.
 *
 * Has the same functionality as AtomCollection, but includes additional functionalities based on periodic boundaries
 */
class PeriodicSystem {
  using masmData =
      std::tuple<AtomCollection, BondOrderCollection, std::unordered_set<unsigned>, std::unordered_map<unsigned, unsigned>>;

 public:
  PeriodicSystem(const PeriodicSystem& other);
  /**
   * @brief Construct a new PeriodicSystem object.
   *
   * All atoms will be created and as type: None and be placed at (0,0,0).
   *
   * @param pbc The periodic boundary conditions.
   * @param N The number of atoms.
   * @param solidStateAtomIndices The indices of atoms that are part of a solid state system.
   */
  explicit PeriodicSystem(const PeriodicBoundaries& pbc = PeriodicBoundaries(), int N = 0,
                          std::unordered_set<unsigned> solidStateAtomIndices = {});
  /**
   * @brief Construct a new PeriodicSystem object
   *
   * @param pbc The periodic boundary conditions.
   * @param elements The elements of the atoms.
   * @param positions The positions of the atoms.
   * @param solidStateAtomIndices The indices of atoms that are part of a solid state system.
   */
  explicit PeriodicSystem(const PeriodicBoundaries& pbc, const ElementTypeCollection& elements,
                          const PositionCollection& positions, std::unordered_set<unsigned> solidStateAtomIndices = {});
  /**
   * @brief Construct a new PeriodicSystem object
   *
   * @param pbc The periodic boundary conditions.
   * @param atoms The atoms.
   * @param solidStateAtomIndices The indices of atoms that are part of a solid state system.
   */
  explicit PeriodicSystem(const PeriodicBoundaries& pbc, AtomCollection atoms,
                          std::unordered_set<unsigned> solidStateAtomIndices = {});
  /**
   * @brief Get an AtomCollection with image atoms that represent the closest image bond partners
   *
   * @param useSolidStateVanDerWaalsBonds Whether bond order collection is generated based Van der Waals radii
   * between solid state atoms instead of nearest neighbors. This will likely overestimate the bonding within a
   * solid state structure, but also avoids that solid state structures may be split up into separate graphs.
   * The solid state image atoms will be checked for actually being in the second shell based on Nearest Neighbors.
   *
   * @return const AtomCollection
   */
  AtomCollection getAtomCollectionWithImages(bool useSolidStateVanDerWaalsBonds = true);
  /**
   * @brief Get the image atoms only.
   *
   * @param useSolidStateVanDerWaalsBonds Whether bond order collection is generated based Van der Waals radii
   * between solid state atoms instead of nearest neighbors. This will likely overestimate the bonding within a
   * solid state structure, but also avoids that solid state structures may be split up into separate graphs.
   * The solid state image atoms will be checked for actually being in the second shell based on Nearest Neighbors.
   *
   * @return AtomCollection
   */
  const AtomCollection& getImageAtoms(bool useSolidStateVanDerWaalsBonds = true);
  /**
   * @brief Get a map pointing from the index of the image atom to the index of the real space atom.
   *
   * @param useSolidStateVanDerWaalsBonds Whether bond order collection is generated based Van der Waals radii
   * between solid state atoms instead of nearest neighbors. This will likely overestimate the bonding within a
   * solid state structure, but also avoids that solid state structures may be split up into separate graphs.
   * The solid state image atoms will be checked for actually being in the second shell based on Nearest Neighbors.
   *
   * @return std::unordered_map<unsigned, unsigned>
   */
  const std::unordered_map<unsigned, unsigned>& getImageAtomsMap(bool useSolidStateVanDerWaalsBonds = true);
  /**
   * @brief Give necessary data for interpret call to ensure valid graph for solid state systems and/or periodic
   * systems.
   *
   * @note all necessary data is constructed if not already present. The bond orders are constructed with the
   * BondDetector.
   *
   * @param useSolidStateVanDerWaalsBonds Whether bond order collection is generated based Van der Waals radii
   * between solid state atoms instead of nearest neighbors. This will likely overestimate the bonding within a
   * solid state structure, but also avoids that solid state structures may be split up into separate graphs.
   * The solid state image atoms will be checked for actually being in the second shell based on Nearest Neighbors.
   *
   * @return masmData A tuple of an AtomCollection including periodic images, a BondOrderCollection including the
   * images, the ignored solidState atom indices and a map specifying the image atom indices
   */
  masmData getDataForMolassemblerInterpretation(bool useSolidStateVanDerWaalsBonds = true);
  /**
   * @brief Give necessary data for interpret call to ensure valid graph for solid state systems and/or periodic
   * systems.
   *
   * @param bondOrderCollection specify bond orders for the system (without image atoms), bonds across periodic
   * boundaries have to be negative.
   * @param removeSolidSecondShell Whether we should remove a second shell of image atoms for solid state atoms
   * based on a nearest neighbor criterion. Recommend if bond orders were generated based on VdW radii.
   *
   * @return masmData A tuple of an AtomCollection including periodic images, a BondOrderCollection including the
   * images, the ignored solidState atom indices and a map specifying the image atom indices
   */
  masmData getDataForMolassemblerInterpretation(const BondOrderCollection& bondOrderCollection,
                                                bool removeSolidSecondShell = false);
  /**
   * @brief Generate image atoms based on BondOrderCollection by BondDetector
   *
   * @note If the AtomCollection of the PeriodicSystem includes (nearly) overlapping atoms, this will most likely fail.
   *
   * @param useSolidStateVanDerWaalsBonds Whether bond order collection is generated based Van der Waals radii
   * between solid state atoms instead of nearest neighbors. This will likely overestimate the bonding within a
   * solid state structure, but also avoids that solid state structures may be split up into separate graphs.
   * The solid state image atoms will be checked for actually being in the second shell based on Nearest Neighbors.
   *
   */
  void constructImageAtoms(bool useSolidStateVanDerWaalsBonds = true);
  /**
   * @brief Generate image atoms based on given BondOrderCollection, negative bond orders indicate bonds through
   * PeriodicBoundaries
   *
   * @note If the AtomCollection of the PeriodicSystem includes (nearly) overlapping atoms, this will most likely fail.
   *
   * @param bondOrders The BondOrderCollection for the PeriodicSystem
   * @param removeSolidSecondShell Whether we should remove a second shell of image atoms for solid state atoms
   * based on a nearest neighbor criterion. Recommend if bond orders were generated based on VdW radii.
   */
  void constructImageAtoms(const BondOrderCollection& bondOrders, bool removeSolidSecondShell = false);
  /**
   * @brief Constructs bond orders of the given AtomCollection with negative bond orders across periodic boundaries
   * based on the BondDetector
   *
   * @param atomCollection The AtomCollection you want the bond orders for.
   * @param periodic Whether periodic boundary conditions should be considered for the bond orders. If they are bonds
   * across periodic boundaries receive a negative bond order.
   * @param useSolidStateVanDerWaalsBonds Whether bond order collection is generated based Van der Waals radii
   * between solid state atoms instead of nearest neighbors. This will likely overestimate the bonding within a
   * solid state structure, but also avoids that solid state structures may be split up into separate graphs.
   * The solid state image atoms will be checked for actually being in the second shell based on Nearest Neighbors.
   *
   * @return BondOrderCollection
   */
  BondOrderCollection constructBondOrders(bool periodic = true, bool useSolidStateVanDerWaalsBonds = true) const;
  /**
   * @brief Takes bond orders, turns them to absolute values and sets all bond orders to negative values, if the bond
   * is spanning across periodic boundaries in the current state of the PeriodicSystem
   *
   * @param bondOrderCollection The template BondOrderCollection specifying the bonds (can be negative or positive).
   */
  void makeBondOrdersAcrossBoundariesNegative(BondOrderCollection& bondOrderCollection) const;
  /**
   * @brief Constructs bond orders of the underlying atoms plus their image atoms based on the BondDetector
   *
   * @note Image atoms are generated on the fly, if they are still missing.
   * @param useSolidStateVanDerWaalsBonds Whether bond order collection is generated based Van der Waals radii
   * between solid state atoms instead of nearest neighbors. This will likely overestimate the bonding within a
   * solid state structure, but also avoids that solid state structures may be split up into separate graphs.
   * The solid state image atoms will be checked for actually being in the second shell based on Nearest Neighbors.
   *
   */
  void constructBondOrdersForImages(bool useSolidStateVanDerWaalsBonds = true);
  /**
   * @brief Constructs bond orders of the underlying atoms plus their image atoms based on the given bond orders.
   *
   * @note In the constructed BondOrderCollection the atoms are bonded both to the real counter part and its closer
   * image.
   *
   * @param bondOrdersWithoutImages A BondOrderCollection of the underlying AtomCollection with negative bond orders
   * signalling bonds across periodic boundaries.
   *
   * @note Image atoms are generated on the fly with the given bond orders, if they are still missing.
   */
  void constructBondOrdersForImages(const BondOrderCollection& bondOrdersWithoutImages);
  /**
   * @brief Get a BondOrderCollection with image atoms that represent the closest image bond partners
   *
   * @note atoms are bonded both to the real counter part and its closer image
   *
   * @param useSolidStateVanDerWaalsBonds Whether bond order collection is generated based Van der Waals radii
   * between solid state atoms instead of nearest neighbors. This will likely overestimate the bonding within a
   * solid state structure, but also avoids that solid state structures may be split up into separate graphs.
   * The solid state image atoms will be checked for actually being in the second shell based on Nearest Neighbors.
   *
   */
  const BondOrderCollection& getBondOrderCollectionWithImages(bool useSolidStateVanDerWaalsBonds = true);
  /**
   * @brief First translates center of mass into the center of the cell and then projects all overhanging atoms
   * back into the cell
   *
   * @note deletes image atoms
   */
  void centerAndTranslateAtomsIntoCell();
  /**
   * @brief Translates all positions of the atoms into the periodic boundaries
   *
   * @note deletes image atoms
   */
  void translateAtomsIntoCell();
  /**
   * @brief Clears all content
   */
  void clear();
  /**
   * @brief Get a PeriodicSystem that represents the primitive cell of the system.
   * @param epsilon The symmetry and equal position precision.
   * @param solidStateOnly Whether you want only a primitive cell of all the solid state atoms, by default false.
   * @return PeriodicSystem The system containing only the primitive cell.
   */
  PeriodicSystem getPrimitiveCellSystem(double epsilon, bool solidStateOnly = false) const;
  /**
   * @brief Performs comparison of the periodic boundaries, solid state atom indices, and contained element
   *   types. It checks the positions also out of order and taking periodic images and layer shifts into account.
   *
   * @param other The PeriodicSystem to compare against
   *
   * @note The positions and PeriodicBoundaries are fuzzy-compared with Eigen's isApprox function and
   *   must therefore not be exactly equal.
   *
   * @return Whether both PeriodicSystems contain the same information
   */
  bool operator==(const PeriodicSystem& other) const;
  //! Negates @see operator ==
  bool operator!=(const PeriodicSystem& other) const;
  /**
   * same logic as @see operator ==
   * allows to set the required accuracy for the fuzzy comparisons
   */
  bool isApprox(const PeriodicSystem& other, double eps = 1e-6) const;
  /**
   * @brief Operator overload, make supersystem
   * @param scalingFactors The scaling factors in x, y, and z
   * @return The super system.
   */
  PeriodicSystem operator*(const Eigen::Vector3i& scalingFactors) const;
  /**
   * @brief Operator overload, make supersystem in place
   * @param scalingFactors The scaling factors in x, y, and z
   * @return The super system.
   */
  PeriodicSystem& operator*=(const Eigen::Vector3i& scalingFactors);
  /**
   * @brief Operator overload, make supersystem
   * @param scalingFactors The scaling factors in x, y, and z
   * @return The super system.
   */
  PeriodicSystem operator*(const std::vector<int>& scalingFactors) const;
  /**
   * @brief Operator overload, make supersystem in place
   * @param scalingFactors The scaling factors in x, y, and z
   * @return The super system.
   */
  PeriodicSystem& operator*=(const std::vector<int>& scalingFactors);
  /**
   * @brief Operator overload, make supersystem
   * @param scalingFactor The scaling factor in all dimensions
   * @return The super system.
   */
  PeriodicSystem operator*(int scalingFactor) const;
  /**
   * @brief Operator overload, make supersystem in place
   * @param scalingFactor The scaling factors in all dimensions
   * @return The super system.
   */
  PeriodicSystem& operator*=(int scalingFactor);

  PeriodicBoundaries pbc;
  AtomCollection atoms;
  std::unordered_set<unsigned> solidStateAtomIndices;

 private:
  bool _lastConstructionUsedVdwBonds = false;
  std::shared_ptr<AtomCollection> _imageAtoms;
  std::shared_ptr<BondOrderCollection> _imageBondOrderCollection;
  // From image atom to base atom index
  // This is necessary because one base atom can have multiple image atoms, but one image atom cannot have
  // multiple real space atoms
  std::unordered_map<unsigned, unsigned> _imageAtomMap;
  // The atoms for which the existing image atoms have been constructed for. Ensure reconstruction if public atom
  // member has been changed in the meantime of last image construction
  AtomCollection _lastImageConstructedAtoms;

  void canonicalize();

  void indicesCheck() const;

  inline void clearImageAtoms() {
    _imageAtoms = nullptr;
    _imageBondOrderCollection = nullptr;
    _imageAtomMap.clear();
  };

  inline bool atomsHaveChangedSinceLastImageConstruction() const {
    return _lastImageConstructedAtoms != atoms;
  }

  void addPotentialImage(int index, const Position& position, bool checkForSecondShell);
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_PERIODICSYSTEM_H_
