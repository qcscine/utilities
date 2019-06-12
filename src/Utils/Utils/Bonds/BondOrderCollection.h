/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_BONDORDERCOLLECTION_H_
#define UTILS_BONDORDERCOLLECTION_H_

#include <Eigen/SparseCore>
#include <cassert>

namespace Scine {
namespace Utils {

/**
 * @class BondOrderCollection BondOrderCollection.h
 * @brief Class defining the bond orders between atoms of some molecular system.
 */
class BondOrderCollection {
 public:
  /// @brief Default Constructor
  BondOrderCollection() = default;
  /**
   * @brief Create an empty collection.
   * @param numberAtoms The number of atoms.
   */
  BondOrderCollection(int numberAtoms) {
    bondOrderMatrix_.resize(numberAtoms, numberAtoms);
  }

  /**
   * @brief Whether this collection is empty, i.e. whether there are any
   *   non-zero bond orders in this collection
   */
  bool empty() const {
    return bondOrderMatrix_.nonZeros() == 0;
  }

  /**
   * @brief Get the Matrix.
   * @return const Eigen::SparseMatrix<double>& The bond order matrix. This is a symmetric
   *   matrix.
   */
  const Eigen::SparseMatrix<double>& getMatrix() const {
    return bondOrderMatrix_;
  }
  /**
   * @brief Set the Matrix.
   * @tparam T Any matrix type.
   * @param bondOrderMatrix The new matrix. This must be a symmetric matrix.
   */
  template<typename T>
  void setMatrix(T&& bondOrderMatrix) {
    assert(bondOrderMatrix.cols() == bondOrderMatrix.rows());
    bondOrderMatrix_ = std::forward<T>(bondOrderMatrix);
  }
  /**
   * @brief Resize the matrix.
   *
   * This operation does not preserve any existing values.
   *
   * @param numberAtoms The new number of atoms.
   */
  void resize(int numberAtoms) {
    bondOrderMatrix_.resize(numberAtoms, numberAtoms);
  }
  /// @brief Set all bond orders to zero and frees the memory.
  void setZero() {
    bondOrderMatrix_.setZero();
    bondOrderMatrix_.data().squeeze();
  }
  /**
   * @brief Get the system size.
   * @tparam Index The type of the index desired. Defaults to integer
   * @return int Returns the number of atoms considered in the BondOrderCollection.
   */
  template<typename Index = int>
  Index getSystemSize() const {
    auto size = bondOrderMatrix_.cols();
    assert(std::numeric_limits<Index>::max() >= size && "Requested matrix size type cannot fit system size");
    return static_cast<Index>(size);
  }
  /**
   * @brief Set the bond order
   *
   * Sets both \f$M_{ij}\f$ and \f$M_{ji}\f$
   *
   * @tparam Index
   * @param i The index of the first atom.
   * @param j The index of the second atom.
   * @param double The new bond order.
   */
  template<typename Index>
  void setOrder(Index i, Index j, double order) {
    bondOrderMatrix_.coeffRef(i, j) = order;
    bondOrderMatrix_.coeffRef(j, i) = order;
  }
  /**
   * @brief Get the Order object
   * @tparam Index
   * @param i The index of the first atom.
   * @param j The index of the second atom.
   * @return double The bond order.
   */
  template<typename Index>
  double getOrder(Index i, Index j) const {
    return bondOrderMatrix_.coeff(i, j);
  }

  /**
   * @brief Checks whether two bond order collections are approximately equal
   *
   * @param other The other bond order collection to compare against
   *
   * @note Uses Eigen's isApprox to fuzzy-compare. Elements in both matrices
   *   need not be exactly the same.
   *
   * @return Whether the matrices are approximately the same
   */
  bool operator==(const BondOrderCollection& other) const {
    return bondOrderMatrix_.isApprox(other.bondOrderMatrix_);
  }

  //! Negates @see operator ==
  bool operator!=(const BondOrderCollection& other) const {
    return !(*this == other);
  }

 private:
  Eigen::SparseMatrix<double> bondOrderMatrix_;
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_BONDORDERCOLLECTION_H_
