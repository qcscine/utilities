/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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

  pybind11::class_<SolventPlacementSettings>(solvation_submodule, "placement_settings")
      .def(pybind11::init<>())
      .def("__repr__", [](const SolventPlacementSettings&) { return "solvent_placement_settings"; })
      .def_readwrite("resolution", &SolventPlacementSettings::resolution)
      .def_readwrite("solvent_offset", &SolventPlacementSettings::solventOffset)
      .def_readwrite("max_distance", &SolventPlacementSettings::maxDistance)
      .def_readwrite("step_size", &SolventPlacementSettings::stepSize)
      .def_readwrite("num_rotamers", &SolventPlacementSettings::numRotamers)
      .def_readwrite("strategic_solvation", &SolventPlacementSettings::strategicSolv)
      .def_readwrite("coverage_threshold", &SolventPlacementSettings::coverageThreshold);

  solvation_submodule.def(
      "solvate",
      pybind11::overload_cast<const AtomCollection&, int, const AtomCollection&, int, int, SolventPlacementSettings>(&solvate),
      pybind11::arg("solute_complex"), pybind11::arg("solute_size"), pybind11::arg("solvent"),
      pybind11::arg("number_solvents"), pybind11::arg("seed"),
      pybind11::arg("placement_settings") = SolventPlacementSettings(),
      "Add systematically a number of one type of "
      "solvent to solute.");

  solvation_submodule.def("solvate_mix", &solvateMix, pybind11::arg("solute_complex"), pybind11::arg("solute_size"),
                          pybind11::arg("solvents"), pybind11::arg("solvent_ratios"), pybind11::arg("number_solvents"),
                          pybind11::arg("seed"), pybind11::arg("placement_settings") = SolventPlacementSettings(),
                          "Add systematically a"
                          "number of different solvents to solute.");

  solvation_submodule.def("solvate_shells", &solvateShells, pybind11::arg("solute_complex"),
                          pybind11::arg("solute_size"), pybind11::arg("solvent"), pybind11::arg("number_shells"),
                          pybind11::arg("seed"), pybind11::arg("placement_settings") = SolventPlacementSettings(),
                          "Add number of one type of solvent in shells to solute.");

  solvation_submodule.def("solvate_shells_mix", &solvateShellsMix, pybind11::arg("solute_complex"),
                          pybind11::arg("solute_size"), pybind11::arg("solvents"), pybind11::arg("solvent_ratios"),
                          pybind11::arg("number_shells"), pybind11::arg("seed"),
                          pybind11::arg("placement_settings") = SolventPlacementSettings(),
                          "Add different"
                          "solvents in shells to solute.");

  solvation_submodule.def("give_solvent_shell_vector", &giveSolventShellVector, pybind11::arg("complex"),
                          pybind11::arg("solute_size"), pybind11::arg("solvent_size_vector"), pybind11::arg("resolution"),
                          pybind11::arg("logger"), pybind11::arg("strategic_solvation") = true,
                          pybind11::arg("threshold") = 1.0, "Analyze a complex and return its solvent shell vector.");

  solvation_submodule.def(
      "transfer_solvent_shell_vector", &transferSolventShellVector, pybind11::arg("shell_vector"),
      "Translate solvent shell vector into one vector containing the size of the solvents in order.");

  solvation_submodule.def("merge_atom_collection_vector", &mergeAtomCollectionVector,
                          pybind11::arg("atom_collection_vector"), "Merges list of atom collections to one atom collection.");

  solvation_submodule.def(
      "merge_solvent_shell_vector", &mergeSolventShellVector, pybind11::arg("shell_vector"),
      "Merge a vector of a vector of atom collections (solvent shell vector) to one atom collection.");

  solvation_submodule.def("merge_solvent_shell_indices_vector", &mergeSolventShellIndices,
                          pybind11::arg("shell_indices_vector"),
                          "Merge a vector of a vector of indices (solvent shell indices vector) to a one flat list.");

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
