/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_BONDDETECTOR_H_
#define UTILS_BONDDETECTOR_H_

#include "Utils/DataStructures/PeriodicBoundaries.h"
#include "Utils/Geometry/PeriodicSystem.h"
#include "Utils/Geometry/Utilities/Distances.h"
#include "Utils/Typenames.h"

namespace Scine {
namespace Utils {

enum class ElementType : unsigned;
class AtomCollection;
class BondOrderCollection;

/**
 * @class BondDetector BondDetector.h
 * @brief Detecting bonds from 3D structures based on interatomic distances using covalent radii.
 *
 * A bond is detected if the distance between two atoms is smaller than the sum of their covalent radii plus 0.4
 * Angstrom. A binary decision on whether a bond exists (resulting in a bond order of 1.0) or not (yielding a bond order
 * of 0.0) is made.
 *
 * The covalent radii were extracted from the Cambridge Structural Database and are available in
 * Utils::BondDetectorRadii.
 *
 * References:\n
 * E. C. Meng, R. A. Lewis, Comput. Chem. 1991, 12, 891-898. [DOI: 10.1002/jcc.540120716]\n
 * C. R. Groom, I. J. Bruno, M. P. Lightfoot and S. C. Ward, Acta Cryst. (2016). B72, 171-179.
 * [DOI: 10.1107/S2052520616003954]\n
 */
class BondDetector {
 public:
  /// @brief Static functions only.
  BondDetector() = delete;
  /**
   * @brief Detect all bonds in an AtomCollection.
   * @param atoms The collection of atoms.
   * @param vanDerWaalsBond Whether van der Waals radii instead of covalent radii should be used.
   * @return BondOrderCollection The bonds.
   */
  static BondOrderCollection detectBonds(const AtomCollection& atoms, bool vanDerWaalsBond = false);
  /**
   * @brief Detect all bonds based on elements and positions.
   * @param elements The collection of elements.
   * @param positions  The collection of positions.
   * @param vanDerWaalsBond Whether van der Waals radii instead of covalent radii should be used.
   * @return BondOrderCollection The bonds.
   */
  static BondOrderCollection detectBonds(const ElementTypeCollection& elements, const PositionCollection& positions,
                                         bool vanDerWaalsBond = false);
  /**
   * @brief Detect all bonds in a PeriodicSystem.
   * @param periodicSystem The Periodic System.
   * @param bondsAcrossBoundariesNegative Whether bonds that are formed across periodic boundaries shall be signalled
   * via negative bond orders. Note that this is very much specific to the given positions.
   * @param vanDerWaalsBond Whether van der Waals radii instead of covalent radii should be used.
   * @return BondOrderCollection The bonds.
   */
  static BondOrderCollection detectBonds(const PeriodicSystem& periodicSystem,
                                         bool bondsAcrossBoundariesNegative = false, bool vanDerWaalsBond = false);
  /**
   * @brief Detect all bonds in an AtomCollection with periodic boundary conditions.
   * @param elements The collection of elements.
   * @param positions  The collection of positions.
   * @param pbc  The Periodic Boundaries.
   * @param bondsAcrossBoundariesNegative Whether bonds that are formed across periodic boundaries shall be signalled
   * via negative bond orders. Note that this is very much specific to the given positions.
   * @param vanDerWaalsBond Whether van der Waals radii instead of covalent radii should be used.
   * @return BondOrderCollection The bonds.
   */
  static BondOrderCollection detectBonds(const AtomCollection& atoms, const PeriodicBoundaries& pbc,
                                         bool bondsAcrossBoundariesNegative = false, bool vanDerWaalsBond = false);
  /**
   * @brief Detect all bonds based on elements and positions with periodic boundary conditions.
   * @param elements The collection of elements.
   * @param positions  The collection of positions.
   * @param pbc  The Periodic Boundaries.
   * @param bondsAcrossBoundariesNegative Whether bonds that are formed across periodic boundaries shall be signalled
   * via negative bond orders. Note that this is very much specific to the given positions.
   * @param vanDerWaalsBond Whether van der Waals radii instead of covalent radii should be used.
   * @return BondOrderCollection The bonds.
   */
  static BondOrderCollection detectBonds(const ElementTypeCollection& elements, const PositionCollection& positions,
                                         const PeriodicBoundaries& pbc, bool bondsAcrossBoundariesNegative = false,
                                         bool vanDerWaalsBond = false);
  /**
   * @brief Check if a bond exists between two atoms.
   * @param e1 The element of the first atom.
   * @param e2 The element of the second atom.
   * @param p1 The position of the first atom.
   * @param p2 The position of the second atom.
   * @param bondsAcrossBoundariesNegative Whether bonds that are formed across periodic boundaries shall be signalled
   * via negative bond orders. Note that this is very much specific to the given positions.
   * @param vanDerWaalsBond Whether van der Waals radii instead of covalent radii should be used.
   * @return true If a bond exists.
   * @return false If no bond exists.
   */
  static bool bondExists(const ElementType& e1, const ElementType& e2, const Position& p1, const Position& p2,
                         bool vanDerWaalsBond = false);
  /**
   * @brief Check if a bond exists between two atoms with periodic boundary conditions.
   * @param e1 The element of the first atom.
   * @param e2 The element of the second atom.
   * @param p1 The position of the first atom.
   * @param p2 The position of the second atom.
   * @param pbc The Periodic Boundaries.
   * @param vanDerWaalsBond Whether van der Waals radii instead of covalent radii should be used.
   * @return true If a bond exists.
   * @return false If no bond exists.
   */
  static bool bondExists(const ElementType& e1, const ElementType& e2, const Position& p1, const Position& p2,
                         const PeriodicBoundaries& pbc, bool vanDerWaalsBond = false);

 private:
  /**
   * @brief Get the covalent radius for a given element.
   * @param e The element.
   * @return double The covalent radius.
   */
  static double getCovalentRadius(ElementType e);
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_BONDDETECTOR_H_
