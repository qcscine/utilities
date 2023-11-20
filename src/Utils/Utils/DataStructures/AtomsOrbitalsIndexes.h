/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_ATOMSORBITALSINDEXES_H
#define UTILS_ATOMSORBITALSINDEXES_H

#include <vector>

namespace Scine {
namespace Utils {

/*!
 * Structure containing the information about AO indexes and their corresponding atom indexes.
 * TODO: separate state from construction (addAtom). In principle nextAtom_ and nextAO_ shouldn't be members of this
 * class.
 */
class AtomsOrbitalsIndexes {
 public:
  explicit AtomsOrbitalsIndexes(int nAtoms = 0);

  /*! Add an atom with nAOs atomic orbitals (for setup). */
  void addAtom(int nAOs);
  /*! Clear the structure. */
  void clear();
  /*! Sets the new number of atoms and resizes members. It is advised to call clear() before. */
  void setSize(int nAtoms);
  /*! Get the total number of atoms. */
  inline int getNAtoms() const {
    return nAtoms_;
  }
  /*! Get the total number of atomic orbitals. */
  inline int getNAtomicOrbitals() const {
    return nAtomicOrbitals_;
  }
  /*! Get the number of atomic orbitals for a given atom index. */
  inline int getNOrbitals(int atomicIndex) const {
    return numberOrbitals_.at(atomicIndex);
  }
  /*! Get the index of the first atomic orbital for a given atom index. */
  inline int getFirstOrbitalIndex(int atomicIndex) const {
    return firstAOIndexes_.at(atomicIndex);
  }

 private:
  int nAtoms_, nAtomicOrbitals_{0};
  std::vector<int> firstAOIndexes_;
  std::vector<int> numberOrbitals_;
  int nextAtom_{0}, nextAO_{0};
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_ATOMSORBITALSINDEXES_H
