/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_ao_to_atom_mapping(pybind11::module& m) {
  pybind11::class_<AtomsOrbitalsIndexes> ao_to_atom_mapping(m, "AOtoAtomMapping");
  ao_to_atom_mapping.def(pybind11::init<int>(), pybind11::arg("nAtoms") = 0, "Initialize a particular number of empty atoms");

  ao_to_atom_mapping.def(
      "get_n_atoms", [](const AtomsOrbitalsIndexes& mapping) -> int { return mapping.getNAtoms(); },
      "Returns the number of atoms.");
  ao_to_atom_mapping.def(
      "get_n_atomic_orbitals", [](const AtomsOrbitalsIndexes& mapping) -> int { return mapping.getNAtomicOrbitals(); },
      "Returns the sum of AOs for all atoms.");
  ao_to_atom_mapping.def(
      "get_n_orbitals",
      [](const AtomsOrbitalsIndexes& mapping, int atomicIndex) -> int { return mapping.getNOrbitals(atomicIndex); },
      "Returns the number of AOs of the atom with the given index.");
  ao_to_atom_mapping.def(
      "get_first_orbital_index",
      [](const AtomsOrbitalsIndexes& mapping, int atomicIndex) -> int {
        return mapping.getFirstOrbitalIndex(atomicIndex);
      },
      "Returns the index of the first AO of the atom with the given index.");
}
