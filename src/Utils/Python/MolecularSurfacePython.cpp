/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Solvation/MolecularSurface.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;
using namespace Scine::Utils::MolecularSurface;

void init_molecular_surface(pybind11::module& m) {
  auto molecular_surface_submodule = m.def_submodule("MolecularSurface");

  pybind11::class_<SurfaceSite> surface_site(molecular_surface_submodule, "SurfaceSite");

  surface_site.def_property_readonly("position", [](const SurfaceSite& o) { return o.position; });
  surface_site.def_property_readonly("normal", [](const SurfaceSite& o) { return o.normal; });

  molecular_surface_submodule.def("get_unpruned_atom_surface", &getUnprunedAtomSurface, pybind11::arg("atom"),
                                  pybind11::arg("atom_surface_points"), "Build unpruned atom surface around an atom");

  molecular_surface_submodule.def("get_pruned_atom_surface", &getPrunedAtomSurface, pybind11::arg("atom_index"),
                                  pybind11::arg("atoms"), pybind11::arg("atom_surface_points"),
                                  "Prune the surface points of one atom in a molecule");

  molecular_surface_submodule.def("get_pruned_molecular_surface", &getPrunedMolecularSurface, pybind11::arg("atoms"),
                                  pybind11::arg("atom_surface_points"), "Gives pruned molecular surface of given molecule.");

  molecular_surface_submodule.def("ray_misses_sphere", &rayMissesSphere, pybind11::arg("surface_site"),
                                  pybind11::arg("sphere_origin"), pybind11::arg("sphere_radius"),
                                  "Checks if surface point misses given sphere");

  molecular_surface_submodule.def("get_visible_molecular_surface", &getVisibleMolecularSurface, pybind11::arg("molecule"),
                                  pybind11::arg("start_index"), pybind11::arg("end_index"), pybind11::arg("resolution") = 64,
                                  "Finds surface sites of molecule in complex which do not 'see' other atoms.");

  molecular_surface_submodule.def(
      "write_surface", &writeSurface, "Write molecular surface into an xyz file with the surface points represented by H atoms.");
}