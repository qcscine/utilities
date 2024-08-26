/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Bonds/BondOrderCollection.h"

namespace Scine {
namespace Utils {

BondOrderCollection::BondOrderCollection(int numberAtoms) {
  bondOrderMatrix_.resize(numberAtoms, numberAtoms);
}

bool BondOrderCollection::empty() const {
  return bondOrderMatrix_.nonZeros() == 0;
}

const Eigen::SparseMatrix<double>& BondOrderCollection::getMatrix() const {
  return bondOrderMatrix_;
}

void BondOrderCollection::resize(int numberAtoms) {
  bondOrderMatrix_.resize(numberAtoms, numberAtoms);
}

void BondOrderCollection::setZero() {
  bondOrderMatrix_.setZero();
  bondOrderMatrix_.data().squeeze();
}

void BondOrderCollection::setToAbsoluteValues() {
  bondOrderMatrix_ = bondOrderMatrix_.cwiseAbs();
}

bool BondOrderCollection::operator==(const BondOrderCollection& other) const {
  return bondOrderMatrix_.isApprox(other.bondOrderMatrix_);
}

bool BondOrderCollection::operator!=(const BondOrderCollection& other) const {
  return !(*this == other);
}

void BondOrderCollection::removeAtomsByIndices(const std::vector<unsigned int>& atomIndices) {
  // Map old to new indices
  std::unordered_map<unsigned int, unsigned int> indexMap;
  unsigned int iNewAtom = 0;
  for (unsigned int iAtom = 0; iAtom < this->bondOrderMatrix_.cols(); ++iAtom) {
    if (std::find(atomIndices.begin(), atomIndices.end(), iAtom) != atomIndices.end()) {
      continue;
    }
    indexMap[iAtom] = iNewAtom;
    ++iNewAtom;
  }
  // Get the bond order values and reset the sparse matrix contents.
  std::vector<Eigen::Triplet<double>> triplets;

  for (const auto& oldToNewI : indexMap) {
    for (Eigen::SparseMatrix<double>::InnerIterator itRow(this->bondOrderMatrix_, oldToNewI.first); itRow; ++itRow) {
      auto oldToNewJ = indexMap.find(itRow.row());
      if (oldToNewJ != indexMap.end()) {
        double value = itRow.value();
        triplets.emplace_back(oldToNewI.second, oldToNewJ->second, value);
      }
    }
  }
  this->bondOrderMatrix_.resize(iNewAtom, iNewAtom);
  this->bondOrderMatrix_.setFromTriplets(triplets.begin(), triplets.end());
}

} /* namespace Utils */
} /* namespace Scine */