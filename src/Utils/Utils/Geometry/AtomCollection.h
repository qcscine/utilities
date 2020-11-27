/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ATOMCOLLECTION_H_
#define UTILS_ATOMCOLLECTION_H_

#include "Utils/Geometry/Atom.h"
#include "Utils/Typenames.h"

namespace Scine {
namespace Utils {

/**
 * @class AtomCollection AtomCollection.h
 * @brief A Collection of Atoms.
 *
 * Has the same functionality as a std::vector<Atom>, but is implemented as a composition of a
 * ElementTypeCollection and a PositionCollection.
 */
class AtomCollection {
 public:
  /**
   * @class AtomCollectionIterator
   * @brief Custom iterator typename for AtomCollection.
   */
  class AtomCollectionIterator : public std::iterator<std::bidirectional_iterator_tag, Atom> {
   public:
    explicit AtomCollectionIterator(AtomCollection const* ac = nullptr, int num = 0);
    AtomCollectionIterator& operator++();
    AtomCollectionIterator operator++(int);
    AtomCollectionIterator& operator--();
    AtomCollectionIterator operator--(int);
    bool operator==(AtomCollectionIterator other) const;
    bool operator!=(AtomCollectionIterator other) const;
    value_type operator*() const;

   private:
    AtomCollection const* ac_;
    int num_;
  };

  using iterator = AtomCollectionIterator;
  /**
   * @brief Construct a new AtomCollection object.
   *
   * All atoms will be created and as type: None and be placed at (0,0,0).
   *
   * @param N The number of atoms.
   */
  explicit AtomCollection(int N = 0);
  /**
   * @brief Construct a new AtomCollection object
   *
   * @param elements The elements of the atoms.
   * @param positions The positions of the atoms.
   */
  AtomCollection(ElementTypeCollection elements, PositionCollection positions);
  /**
   * @brief Update all Elements.
   *
   * The new data must have the correct size,
   * this function does not resize the data objects.
   *
   * @param elements The new collection of Elements.
   */
  void setElements(ElementTypeCollection elements);
  /**
   * @brief Update all positions
   *
   * The new data must have the correct size,
   * this function does not resize the data objects.
   *
   * @param positions The new positions.
   */
  void setPositions(PositionCollection positions);
  /**
   * @brief Get the Elements object
   *
   * @return const ElementTypeCollection&
   */
  const ElementTypeCollection& getElements() const;
  /**
   * @brief Get the Positions object
   *
   * @return const PositionCollection&
   */
  const PositionCollection& getPositions() const;
  /**
   * @brief Set the Element of one existing atom.
   *
   * This function can not access/create new atoms that are not present yet.
   *
   * @param i The index of the atom.
   * @param e The new Element.
   */
  void setElement(int i, ElementType e);
  /**
   * @brief Set the Position of one existing atom.
   *
   * This function can not access/create new atoms that are not present yet.
   *
   * @param i The index of the atom.
   * @param p The new Position.
   */
  void setPosition(int i, const Position& p);
  /**
   * @brief Get the Element of a single atom.
   * @param i The index of the atom.
   * @return ElementType The Element of atom i.
   */
  ElementType getElement(int i) const;
  /**
   * @brief Get the Position of a a single atom.
   * @param i The index of the atom.
   * @return Position The Position of atom i.
   */
  Position getPosition(int i) const;
  /**
   * @brief Getter for thw collection size.
   * @return int Returns the number of atoms in the collection.
   */
  int size() const;
  /**
   * @brief The start of the collection iterator.
   * @return iterator The iterator.
   */
  iterator begin() const;
  /**
   * @brief The end of the collection iterator.
   * @return iterator The iterator.
   */
  iterator end() const;
  /**
   * @brief Clears all content and resizes data objects to size 0;
   */
  void clear();
  /**
   * @brief Resizes the data obejects to the given size.
   * @param n The number of atoms to resize to.
   */
  void resize(int n);
  /**
   * @brief Appends one atom to the collectio via copy.
   * @param atom The atom.
   */
  void push_back(const Atom& atom);
  /**
   * @brief Operator overload, getter by position.
   * @param i The position.
   * @return Atom That atom at position i.
   */
  Atom operator[](int i) const;
  /**
   * @brief Getter by position.
   * @param i The position.
   * @return Atom That atom at position i.
   */
  Atom at(int i) const;
  /**
   * @brief Performs in-order comparison of both the contained element
   *   types and positions.
   *
   * @param other The AtomCollection to compare against
   *
   * @note The positions are fuzzy-compared with Eigen's isApprox function and
   *   must therefore not be exactly equal.
   *
   * @return Whether both atomcollection contain the same information
   */
  bool operator==(const AtomCollection& other) const;

  //! Negates @see operator ==
  bool operator!=(const AtomCollection& other) const;
  /**
   * @brief Operator overload, combine two atom collections.
   * @param other The other atom collection, appended to the first one.
   * @return The combined atom collection.
   */
  AtomCollection operator+(const AtomCollection& other) const;
  /**
   * @brief Operator overload, append atom collection other.
   * @param other The atom collection which shall be appended.
   * @return The combined atom collection.
   */
  AtomCollection& operator+=(const AtomCollection& other);

 private:
  ElementTypeCollection elements_;
  PositionCollection positions_;
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_ATOMCOLLECTION_H_
