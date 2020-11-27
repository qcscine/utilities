/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Solvation/SoluteSolventComplex.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;
using namespace Scine::Utils::SoluteSolventComplex;

void init_solute_solvent_complex(pybind11::module& m) {
  auto solvation_submodule = m.def_submodule("solvation");

  solvation_submodule.def(
      "solvate",
      pybind11::overload_cast<const AtomCollection&, int, const AtomCollection&, int, int, int, double, double, double, int, bool>(
          &solvate),
      pybind11::arg("solute_complex"), pybind11::arg("solute_size"), pybind11::arg("solvent"),
      pybind11::arg("number_solvents"), pybind11::arg("seed"), pybind11::arg("resolution") = 32,
      pybind11::arg("solvent_offset") = 0.0, pybind11::arg("max_distance") = 10.0, pybind11::arg("step_size") = 0.25,
      pybind11::arg("number_rotamers") = 3, pybind11::arg("strategic_solvation") = false,
      "Add systematically a number of solvents to solute.");

  solvation_submodule.def("solvate_shells", &solvateShells, pybind11::arg("solute_complex"),
                          pybind11::arg("solute_size"), pybind11::arg("solvent"), pybind11::arg("number_shells"),
                          pybind11::arg("seed"), pybind11::arg("resolution") = 32,
                          pybind11::arg("solvent_offset") = 0.0, pybind11::arg("max_distance") = 10.0,
                          pybind11::arg("step_size") = 0.25, pybind11::arg("number_rotamers") = 3,
                          pybind11::arg("strategic_solvation") = false, "Add number of solvent shells to solute.");

  solvation_submodule.def("merge_atom_collection_vector", &mergeAtomCollectionVector,
                          pybind11::arg("atom_collection_vector"), "Merges list of atom collections to one atom collection.");

  solvation_submodule.def(
      "merge_solvent_shell_vector", &mergeSolventShellVector, pybind11::arg("shell_vector"),
      "Merge a vector of a vector of atom collections (solvent shell vector) to one atom collection.");

  solvation_submodule.def("check_distances", &checkDistances, pybind11::arg("molecule_1"), pybind11::arg("molecule_2"),
                          "Check if two atom collections overlap with their VdW spheres.");

  solvation_submodule.def("solvation_strategy", &solvationStrategy,
                          "Solvation strategy for faster building of solute - solvent complexes.");

  solvation_submodule.def("arrange", &arrange, pybind11::arg("surface_point_1"), pybind11::arg("surface_normal_1"),
                          pybind11::arg("surface_point_2"), pybind11::arg("surface_normal_2"),
                          pybind11::arg("molecule_2"), pybind11::arg("distance"),
                          "Arrange one atom collection such that the two positions given face each other.");

  solvation_submodule.def("add", &add, pybind11::arg("complex"), pybind11::arg("additive"),
                          pybind11::arg("complex_surface_site"), pybind11::arg("additive_surface_site"),
                          pybind11::arg("min_distance"), pybind11::arg("max_distance"),
                          pybind11::arg("increment_distance") = 0.25, pybind11::arg("number_rotation_attempts") = 3,
                          "Add additive to given complex at given surface site of the complex.");
}
