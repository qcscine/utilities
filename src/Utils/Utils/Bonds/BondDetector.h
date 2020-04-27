/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_BONDDETECTOR_H_
#define UTILS_BONDDETECTOR_H_

#include "Utils/Typenames.h"

namespace Scine {
namespace Utils {

enum class ElementType : unsigned;
class AtomCollection;
class BondOrderCollection;

/**
 * @class BondDetector BondDetector.h
 * @brief Detecting bonds from 3D structures using covalent radii.
 *
 * The radii are stored in Scine::Utils::BondDetectorRadii .
 *
 * References:
 * - DOI: 10.1002/jcc.540120716
 * - DOI: 10.1186/1758-2946-4-26
 * - DOI: 10.1002/jcc.24309
 */
class BondDetector {
 public:
  /// @brief Static functions only.
  BondDetector() = delete;
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
  static bool bondExists(const ElementType& e1, const ElementType& e2, const Position& p1, const Position& p2);

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
