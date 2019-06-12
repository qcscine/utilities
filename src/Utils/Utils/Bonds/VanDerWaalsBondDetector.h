/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_VANDERWAALSBONDDETECTOR_H_
#define UTILS_VANDERWAALSBONDDETECTOR_H_

#include "Utils/Typenames.h"

namespace Scine {
namespace Utils {
class AtomCollection;
class BondOrderCollection;
enum class ElementType;

/**
 * @class VanDerWaalsBondDetector VanDerWaalsBondDetector.h
 * @brief Implements the function to find bonds from a set of elements and positions, using the van der Waals radii.
 */
class VanDerWaalsBondDetector {
 public:
  /// @brief Static functions only.
  VanDerWaalsBondDetector() = delete;
  /**
   * @brief Detect all bonds in an AtomCollection.
   * @param atoms The collection of atoms.
   * @return BondOrderCollection The bonds.
   */
  static BondOrderCollection detectBonds(const AtomCollection& atoms);
  /**
   * @brief Detect all bonds based on elements and positions.
   * @param elements The collection of elements.
   * @param positions  The collection of positions.
   * @return BondOrderCollection The bonds.
   */
  static BondOrderCollection detectBonds(const ElementTypeCollection& elements, const PositionCollection& positions);
  /**
   * @brief Check if a bond exists between two atoms.
   * @param e1 The element of the first atom.
   * @param e2 The element of the second atom.
   * @param p1 The position of the first atom.
   * @param p2 The position of the second atom.
   * @return true If a bond exists.
   * @return false If no bond exists.
   */
  static bool vdwBondExists(const ElementType& e1, const ElementType& e2, const Position& p1, const Position& p2);
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_VANDERWAALSBONDDETECTOR_H_
