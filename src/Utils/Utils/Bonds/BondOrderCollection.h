/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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
  /// @brief Set all bond orders to their absolute value.
  void setToAbsoluteValues() {
    bondOrderMatrix_ = bondOrderMatrix_.cwiseAbs();
  }
  /**
   * @brief Get the system size.
   * @tparam Index The type of the index desired. Defaults to integer
   * @return int Returns the number of atoms considered in the BondOrderCollection.
   */
  template<typename Index = int>
  Index getSystemSize() const {
    auto size = bondOrderMatrix_.cols();
#ifndef NDEBUG
    using Common = std::common_type_t<Eigen::Index, Index>;
    assert(static_cast<Common>(std::numeric_limits<Index>::max()) >= static_cast<Common>(size) &&
           "Requested matrix size type cannot fit system size");
#endif
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
    rangeCheck(i, j);
    bondOrderMatrix_.coeffRef(i, j) = order;
    bondOrderMatrix_.coeffRef(j, i) = order;
    if (std::fabs(order) < 1e-12) {
      bondOrderMatrix_.prune(0.0); // gets rid of zero entry to avoid wrong empty() evaluation
    }
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
    rangeCheck(i, j);
    return bondOrderMatrix_.coeff(i, j);
  }
  /**
   * @brief Get bond partner of an index based on a minimum bond threshold
   * @tparam Index
   * @param i The index of whom we want the bond partner.
   * @param threshold The threshold above which a bond exists.
   * @return std::vector<Index> The list of bond partners.
   */
  template<typename Index>
  std::vector<Index> getBondPartners(Index i, double threshold = 0.5) const {
    rangeCheck(i, i);
    std::vector<Index> result;
    for (Index j = 0; j < getSystemSize(); ++j) {
      if (bondOrderMatrix_.coeff(i, j) > threshold) {
        result.push_back(j);
      }
    }
    return result;
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

  template<typename Index>
  inline void rangeCheck(Index i, Index j) const {
    if (i >= getSystemSize<Index>()) {
      throw std::runtime_error("The given index " + std::to_string(i) + " is too big for this BondOrderCollection.");
    }
    if (j >= getSystemSize<Index>()) {
      throw std::runtime_error("The given index " + std::to_string(j) + " is too big for this BondOrderCollection.");
    }
    if (i < 0 || j < 0) {
      throw std::runtime_error(
          "It is not possible to access an element of a BondOrderCollection with a negative index.");
    }
  }
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_BONDORDERCOLLECTION_H_
