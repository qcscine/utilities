/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ATOM_H_
#define UTILS_ATOM_H_

#include "Utils/Geometry/ElementTypes.h"
#include "Utils/Typenames.h"
#include <utility>

namespace Scine {
namespace Utils {

/**
 * @class Atom Atom.h
 * @brief A basic atom in the chemical sense.
 */
class Atom {
 public:
  /**
   * @brief Construct a new Atom:: Atom object
   * @param e The ElementType. [Default: None]
   * @param p The Position. [Default: (0.0,0.0,0.0)]
   */
  explicit Atom(ElementType e = ElementType::none, Position p = Position(0, 0, 0));
  /**
   * @brief Getter for the ElementType.
   * @return ElementType Returns the type of the atom.
   */
  ElementType getElementType() const;
  /**
   * @brief Getter for the Position.
   * @return const Position& Returns the Position of the atom.
   */
  const Position& getPosition() const;
  /**
   * @brief Setter for the ElementType.
   * @param e The ElementType.
   */
  void setElementType(ElementType e);
  /**
   * @brief Setter for the Position.
   * @param pos The Position.
   */
  void setPosition(Position pos);

 private:
  ElementType element_;
  Position position_;
};

inline Atom::Atom(ElementType e, Position p) : element_(e), position_(std::move(p)) {
}

inline ElementType Atom::getElementType() const {
  return element_;
}

inline void Atom::setPosition(Position pos) {
  position_ = std::move(pos);
}

inline void Atom::setElementType(ElementType e) {
  element_ = e;
}

inline const Position& Atom::getPosition() const {
  return position_;
}

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_ATOM_H_
