/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry/AtomCollection.h"

namespace Scine {
namespace Utils {

AtomCollection::AtomCollection(int N) : elements_(N), positions_(N, 3) {
  positions_.setZero();
}

AtomCollection::AtomCollection(ElementTypeCollection elements, PositionCollection positions)
  : elements_(std::move(elements)), positions_(std::move(positions)) {
  assert(elements_.size() == positions_.rows());
}

void AtomCollection::setElements(ElementTypeCollection elements) {
  assert(elements.size() == size());
  elements_ = std::move(elements);
}

void AtomCollection::setPositions(PositionCollection positions) {
  assert(positions.rows() == size());
  positions_ = std::move(positions);
}

void AtomCollection::clear() {
  elements_.clear();
  positions_.resize(0, 3);
}

void AtomCollection::resize(int n) {
  elements_.resize(n);
  positions_.resize(n, 3);
}

AtomCollection::iterator AtomCollection::begin() const {
  return AtomCollection::iterator(this, 0);
}

AtomCollection::iterator AtomCollection::end() const {
  return AtomCollection::iterator(this, size());
}

void AtomCollection::push_back(const Atom& atom) {
  elements_.push_back(atom.getElementType());
  positions_.conservativeResize(positions_.rows() + 1, 3);
  positions_.row(positions_.rows() - 1) = atom.getPosition();
}

bool AtomCollection::operator==(const AtomCollection& other) const {
  return (elements_ == other.elements_ && positions_.isApprox(other.positions_));
}

bool AtomCollection::operator!=(const AtomCollection& other) const {
  return !(*this == other);
}

const ElementTypeCollection& AtomCollection::getElements() const {
  return elements_;
}

const PositionCollection& AtomCollection::getPositions() const {
  return positions_;
}

void AtomCollection::setElement(int i, ElementType e) {
  assert(0 <= i && i < size());
  elements_[i] = e;
}

void AtomCollection::setPosition(int i, const Position& p) {
  assert(0 <= i && i < size());
  positions_.row(i) = p;
}

ElementType AtomCollection::getElement(int i) const {
  assert(0 <= i && i < size());
  return elements_[i];
}

Position AtomCollection::getPosition(int i) const {
  assert(0 <= i && i < size());
  return positions_.row(i);
}

int AtomCollection::size() const {
  return elements_.size();
}

AtomCollection AtomCollection::operator+(const AtomCollection& other) const {
  AtomCollection merge = *this;
  for (const auto& atom : other) {
    merge.push_back(atom);
  }
  return merge;
}

AtomCollection& AtomCollection::operator+=(const AtomCollection& other) {
  AtomCollection merge = *this + other;
  *this = merge;
  return *this;
}

Atom AtomCollection::operator[](int i) const {
  return Atom(elements_[i], positions_.row(i));
}

Atom AtomCollection::at(int i) const {
  return Atom(elements_.at(i), positions_.row(i));
}

AtomCollection::AtomCollectionIterator AtomCollection::AtomCollectionIterator::operator++(int) {
  AtomCollectionIterator retval;
  retval = *this;
  ++(*this);
  return retval;
}

AtomCollection::AtomCollectionIterator& AtomCollection::AtomCollectionIterator::operator++() {
  ++num_;
  return *this;
}

AtomCollection::AtomCollectionIterator AtomCollection::AtomCollectionIterator::operator--(int) {
  AtomCollectionIterator retval = *this;
  --(*this);
  return retval;
}

AtomCollection::AtomCollectionIterator& AtomCollection::AtomCollectionIterator::operator--() {
  --num_;
  return *this;
}

Atom AtomCollection::AtomCollectionIterator::operator*() const {
  return (*ac_)[num_];
}

bool AtomCollection::AtomCollectionIterator::operator==(AtomCollectionIterator other) const {
  return num_ == other.num_;
}

bool AtomCollection::AtomCollectionIterator::operator!=(AtomCollectionIterator other) const {
  return !(*this == other);
}

AtomCollection::AtomCollectionIterator::AtomCollectionIterator(AtomCollection const* ac, int num) : ac_(ac), num_(num) {
}

} /* namespace Utils */
} /* namespace Scine */
