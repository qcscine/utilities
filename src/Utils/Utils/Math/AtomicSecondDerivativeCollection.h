/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_ATOMICSECONDDERIVATIVECOLLECTION_H
#define UTILS_MATH_ATOMICSECONDDERIVATIVECOLLECTION_H

#include "AutomaticDifferentiation/Second3D.h"
#include <Eigen/Core>
#include <algorithm>
#include <list>
#include <vector>

namespace Scine {
namespace Utils {

using SecondDerivative = AutomaticDifferentiation::Second3D;

/**
 * @class AtomicSecondDerivativeCollection AtomicSecondDerivativeCollection.h
 * @brief Container class for second derivatives.
 *
 *      Mainly an override of std::vector<SecondDerivative>.
 *      Generally used for storing the gradients of a molecular structure.
 */
class AtomicSecondDerivativeCollection {
 public:
  using container = std::vector<SecondDerivative>;
  using iterator = container::iterator;
  using const_iterator = container::const_iterator;
  using reference = container::reference;
  using const_reference = container::const_reference;

  /**
   * @brief Constructor for this container class.
   */
  explicit AtomicSecondDerivativeCollection(int N = 0, SecondDerivative d = SecondDerivative());

  iterator begin();
  const_iterator begin() const;
  iterator end();
  const_iterator end() const;

  /**
   * @brief Getter for atomic hessian
   * @param i The index of the atom.
   * @return Eigenmatrix containing the atomic hessian of atom i.
   */
  Eigen::Matrix3d getAtomicHessian(int i) const;
  /**
   * @brief Returns a list of all atomic hessian Eigenmatrices
   * The atomic Hessian of a given atom i within a molecule is the 3x3 sub block of the full
   * Hessian matrix of the molecule which contains only second derivatives with respect to the
   * coordinates of atom i.
   */
  std::list<Eigen::Matrix3d> getAtomicHessians() const;
  /**
   * @brief Clearing the vector that holds the Second3D objects.
   */
  void clear();
  /**
   * @brief Wrapper around standard resizing function for vector that holds the Second3D objects.
   */
  void resize(int n);
  /**
   * @brief Wrapper around push back function for vector of Second3D objects.
   */
  void push_back(const SecondDerivative& d);
  /**
   * @brief Wrapper around push back function for vector of Second3D objects. This time the Second3d object
   *        is given as an rvalue reference.
   */
  void push_back(SecondDerivative&& d);
  reference operator[](int i);
  const_reference operator[](int i) const;
  reference at(int i);
  const_reference at(int i) const;
  iterator erase(iterator position);

  /**
   * @brief Returns true if the vector that holds the Second3D objects is empty.
   */
  bool empty() const;
  /**
   * @brief Returns the size of the vector holding the Second3D objects (dimension of the container).
   */
  int size() const;

  /**
   * @brief Multiplication assignment operator (for example for unit conversion).
   */
  const AtomicSecondDerivativeCollection& operator*=(double f);
  /**
   * @brief Division assignment operator (for example for unit conversion).
   */
  const AtomicSecondDerivativeCollection& operator/=(double f);
  /**
   * @brief Add another AtomicSecondDerivativeCollection.
   */
  const AtomicSecondDerivativeCollection& operator+=(const AtomicSecondDerivativeCollection& dc);
  /**
   * @brief Multiplication operator (for example for unit conversion).
   */
  AtomicSecondDerivativeCollection operator*(double f) const;
  /**
   * @brief Division operator (for example for unit conversion).
   */
  AtomicSecondDerivativeCollection operator/(double f) const;
  /**
   * @brief Reset all the derivatives to zero.
   */
  void setZero();

 private:
  // vector holding the derivatives (container = std::vector<SecondDerivative>)
  container derivativeVector_;
};

// Constructor
inline AtomicSecondDerivativeCollection::AtomicSecondDerivativeCollection(int N, SecondDerivative d)
  : derivativeVector_(static_cast<container::size_type>(N), d) {
}

/*
 *
 * Inline implementation of all the wrapper functions around the std::vector<SecondDerivative>
 *
 */
inline AtomicSecondDerivativeCollection::iterator AtomicSecondDerivativeCollection::begin() {
  return derivativeVector_.begin();
}

inline AtomicSecondDerivativeCollection::const_iterator AtomicSecondDerivativeCollection::begin() const {
  return derivativeVector_.begin();
}

inline AtomicSecondDerivativeCollection::iterator AtomicSecondDerivativeCollection::end() {
  return derivativeVector_.end();
}

inline AtomicSecondDerivativeCollection::const_iterator AtomicSecondDerivativeCollection::end() const {
  return derivativeVector_.end();
}

inline Eigen::Matrix3d AtomicSecondDerivativeCollection::getAtomicHessian(int i) const {
  assert(0 <= i && i < size());
  Eigen::Matrix3d atomichessianMatrix;
  atomichessianMatrix(0, 0) = derivativeVector_[i].XX();
  atomichessianMatrix(1, 0) = derivativeVector_[i].YX();
  atomichessianMatrix(2, 0) = derivativeVector_[i].ZX();
  atomichessianMatrix(0, 1) = derivativeVector_[i].XY();
  atomichessianMatrix(1, 1) = derivativeVector_[i].YY();
  atomichessianMatrix(2, 1) = derivativeVector_[i].ZY();
  atomichessianMatrix(0, 2) = derivativeVector_[i].XZ();
  atomichessianMatrix(1, 2) = derivativeVector_[i].YZ();
  atomichessianMatrix(2, 2) = derivativeVector_[i].ZZ();
  return atomichessianMatrix;
}

inline std::list<Eigen::Matrix3d> AtomicSecondDerivativeCollection::getAtomicHessians() const {
  std::list<Eigen::Matrix3d> atomichessianList;
  for (int i = 0; i < size(); ++i) {
    atomichessianList.push_back(AtomicSecondDerivativeCollection::getAtomicHessian(i));
  }
  return atomichessianList;
}

inline void AtomicSecondDerivativeCollection::clear() {
  derivativeVector_.clear();
}

inline void AtomicSecondDerivativeCollection::resize(int n) {
  derivativeVector_.resize(static_cast<container::size_type>(n));
}

inline void AtomicSecondDerivativeCollection::push_back(const SecondDerivative& d) {
  derivativeVector_.push_back(d);
}

inline void AtomicSecondDerivativeCollection::push_back(SecondDerivative&& d) {
  derivativeVector_.push_back(std::move(d));
}

inline AtomicSecondDerivativeCollection::reference AtomicSecondDerivativeCollection::operator[](int i) {
  return derivativeVector_[i];
}

inline AtomicSecondDerivativeCollection::reference AtomicSecondDerivativeCollection::at(int i) {
  return derivativeVector_.at(static_cast<container::size_type>(i));
}

inline AtomicSecondDerivativeCollection::const_reference AtomicSecondDerivativeCollection::at(int i) const {
  return derivativeVector_.at(static_cast<container::size_type>(i));
}

inline AtomicSecondDerivativeCollection::const_reference AtomicSecondDerivativeCollection::operator[](int i) const {
  return derivativeVector_[i];
}

inline AtomicSecondDerivativeCollection::iterator AtomicSecondDerivativeCollection::erase(iterator position) {
  return derivativeVector_.erase(position);
}

inline bool AtomicSecondDerivativeCollection::empty() const {
  return derivativeVector_.empty();
}

inline int AtomicSecondDerivativeCollection::size() const {
  return static_cast<int>(derivativeVector_.size());
}

/*
 *
 * Inline implementation of the operators
 *
 */

inline const AtomicSecondDerivativeCollection& AtomicSecondDerivativeCollection::operator*=(double f) {
  for (auto& d : derivativeVector_) {
    d *= f;
  }
  return *this;
}

inline const AtomicSecondDerivativeCollection& AtomicSecondDerivativeCollection::operator/=(double f) {
  for (auto& d : derivativeVector_) {
    d /= f;
  }
  return *this;
}

inline const AtomicSecondDerivativeCollection& AtomicSecondDerivativeCollection::
operator+=(const AtomicSecondDerivativeCollection& dc) {
  // Check that the two collections are of the same size before they are added.
  assert(size() == dc.size() && "Adding AtomicSecondDerivativeCollections of different sizes.");
  const auto s = size();
  for (int i = 0; i < s; ++i) {
    derivativeVector_[i] += dc[i];
  }
  return *this;
}

inline AtomicSecondDerivativeCollection AtomicSecondDerivativeCollection::operator*(double f) const {
  AtomicSecondDerivativeCollection d(size());
  for (int i = 0; i < size(); ++i) {
    d[i] = derivativeVector_[i] * f;
  }
  return d;
}

inline AtomicSecondDerivativeCollection AtomicSecondDerivativeCollection::operator/(double f) const {
  AtomicSecondDerivativeCollection d(size());
  for (int i = 0; i < size(); ++i) {
    d[i] = derivativeVector_[i] / f;
  }
  return d;
}

// Set all derivatives to Zero.
inline void AtomicSecondDerivativeCollection::setZero() {
  for_each(derivativeVector_.begin(), derivativeVector_.end(), [](SecondDerivative& g) { g = SecondDerivative(); });
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_ATOMICSECONDDERIVATIVECOLLECTION_H
