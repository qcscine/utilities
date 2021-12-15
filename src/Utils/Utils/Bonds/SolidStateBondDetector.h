/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_SOLIDSTATEBONDDETECTOR_H_
#define UTILS_SOLIDSTATEBONDDETECTOR_H_

#include "BondDetector.h"
#include "Utils/DataStructures/PeriodicBoundaries.h"
#include "Utils/Geometry/Utilities/Distances.h"
#include "Utils/Typenames.h"
#include <unordered_set>

namespace Scine {
namespace Utils {

enum class ElementType : unsigned;
class AtomCollection;
class BondOrderCollection;

/**
 * @class SolidStateBondDetector SolidStateBondDetector.h
 * @brief Detecting bonds from 3D structures based on interatomic distances using nearest neighbors and covalent radii.
 *
 * This detector can handle both pure solid state structures and heterogeneous systems.
 * For the given solid state atoms, nearest neighbor bond orders are used with a margin of 0.1 Angstrom. If it is a
 * heterogeneous system the conventional BondDetector (see BondDetector.h) is used to determine any bond involving at
 * least one non-solid state atom.
 * A binary decision on whether a bond exists (resulting in a bond order of 1.0) or not (yielding a bond order
 * of 0.0) is made.
 */
class SolidStateBondDetector {
 public:
  /// @brief Static functions only.
  SolidStateBondDetector() = delete;
  /**
   * @brief Detect all bonds in an AtomCollection.
   * @param atoms The collection of atoms.
   * @param solidStateIndices The indices of atoms that are part of the solid phase.
   * @return BondOrderCollection The bonds.
   */
  static BondOrderCollection detectBonds(const AtomCollection& atoms, const std::unordered_set<unsigned>& solidStateIndices);
  /**
   * @brief Detect all bonds based on elements and positions.
   * @param elements The collection of elements.
   * @param positions  The collection of positions.
   * @param solidStateIndices The indices of atoms that are part of the solid phase.
   * @return BondOrderCollection The bonds.
   */
  static BondOrderCollection detectBonds(const ElementTypeCollection& elements, const PositionCollection& positions,
                                         const std::unordered_set<unsigned>& solidStateIndices);
  /**
   * @brief Detect all bonds in a PeriodicSystem.
   * @param periodicSystem The Periodic System, which also holds the solidStateIndices
   * @param bondsAcrossBoundariesNegative Whether bonds that are formed across periodic boundaries shall be signalled
   * via negative bond orders. Note that this is very much specific to the given positions.
   * @return BondOrderCollection The bonds.
   */
  static BondOrderCollection detectBonds(const PeriodicSystem& periodicSystem, bool bondsAcrossBoundariesNegative = false);
  /**
   * @brief Detect all bonds in an AtomCollection with periodic boundary conditions.
   * @param elements The collection of elements.
   * @param positions  The collection of positions.
   * @param pbc  The Periodic Boundaries.
   * @param solidStateIndices The indices of atoms that are part of the solid phase.
   * @param bondsAcrossBoundariesNegative Whether bonds that are formed across periodic boundaries shall be signalled
   * via negative bond orders. Note that this is very much specific to the given positions.
   * @return BondOrderCollection The bonds.
   */
  static BondOrderCollection detectBonds(const AtomCollection& atoms, const PeriodicBoundaries& pbc,
                                         const std::unordered_set<unsigned>& solidStateIndices,
                                         bool bondsAcrossBoundariesNegative = false);
  /**
   * @brief Detect all bonds based on elements and positions with periodic boundary conditions.
   * @param elements The collection of elements.
   * @param positions  The collection of positions.
   * @param pbc  The Periodic Boundaries.
   * @param solidStateIndices The indices of atoms that are part of the solid phase.
   * @param bondsAcrossBoundariesNegative Whether bonds that are formed across periodic boundaries shall be signalled
   * via negative bond orders. Note that this is very much specific to the given positions.
   * @return BondOrderCollection The bonds.
   */
  static BondOrderCollection detectBonds(const ElementTypeCollection& elements, const PositionCollection& positions,
                                         const PeriodicBoundaries& pbc, const std::unordered_set<unsigned>& solidStateIndices,
                                         bool bondsAcrossBoundariesNegative = false);
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_SOLIDSTATEBONDDETECTOR_H_
