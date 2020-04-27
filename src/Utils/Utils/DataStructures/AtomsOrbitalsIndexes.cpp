/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AtomsOrbitalsIndexes.h"

namespace Scine {
namespace Utils {

AtomsOrbitalsIndexes::AtomsOrbitalsIndexes(int nAtoms)
  : nAtoms_(nAtoms), firstAOIndexes_(nAtoms), numberOrbitals_(nAtoms) {
}

void AtomsOrbitalsIndexes::addAtom(int nAOs) {
  firstAOIndexes_[nextAtom_] = nextAO_;
  numberOrbitals_[nextAtom_] = nAOs;
  nAtomicOrbitals_ += nAOs;

  nextAO_ += nAOs;
  ++nextAtom_;
}

void AtomsOrbitalsIndexes::clear() {
  nAtoms_ = 0;
  nAtomicOrbitals_ = 0;
  firstAOIndexes_.clear();
  numberOrbitals_.clear();
  nextAO_ = 0;
  nextAtom_ = 0;
}

void AtomsOrbitalsIndexes::setSize(int nAtoms) {
  nAtoms_ = nAtoms;
  firstAOIndexes_.resize(nAtoms);
  numberOrbitals_.resize(nAtoms);
}

} // namespace Utils
} // namespace Scine
