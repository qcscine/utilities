/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/PeriodicSystem.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_periodic_system(pybind11::module& m) {
  pybind11::class_<PeriodicSystem> periodic_system(m, "PeriodicSystem");
  periodic_system.def(pybind11::init<const PeriodicBoundaries&, int, std::unordered_set<unsigned>>(), pybind11::arg("pbc"),
                      pybind11::arg("N") = 0, pybind11::arg("solid_state_atom_indices") = std::unordered_set<unsigned>(),
                      "Initialize a particular number of empty atoms");
  periodic_system.def(
      pybind11::init<const PeriodicBoundaries&, const ElementTypeCollection&, const PositionCollection&, std::unordered_set<unsigned>>(),
      pybind11::arg("pbc"), pybind11::arg("elements"), pybind11::arg("positions"),
      pybind11::arg("solid_state_atom_indices") = std::unordered_set<unsigned>(),
      "Initialize from element types and positions");
  periodic_system.def(pybind11::init<const PeriodicBoundaries&, AtomCollection, std::unordered_set<unsigned>>(),
                      pybind11::arg("pbc"), pybind11::arg("atoms"),
                      pybind11::arg("solid_state_atom_indices") = std::unordered_set<unsigned>(), "Initialize from atoms");

  // public members
  periodic_system.def_readwrite("pbc", &PeriodicSystem::pbc, "The periodic boundary conditions");
  periodic_system.def_readwrite("atoms", &PeriodicSystem::atoms, "The atoms");
  periodic_system.def_readwrite("solid_state_atom_indices", &PeriodicSystem::solidStateAtomIndices,
                                "The indices of solid state atoms");

  // private members getters
  periodic_system.def("get_atom_collection_with_images", &PeriodicSystem::getAtomCollectionWithImages,
                      "Get the atoms plus the image atoms");
  periodic_system.def("get_image_atoms", &PeriodicSystem::getImageAtoms, "Get the image atoms only");
  periodic_system.def("get_image_atoms_map", &PeriodicSystem::getImageAtomsMap,
                      "Get a map pointing from the index of the image atom to the index of the real space atom");

  // methods surrounding image atoms
  periodic_system.def("get_data_for_molassembler_interpretation",
                      pybind11::overload_cast<>(&PeriodicSystem::getDataForMolassemblerInterpretation),
                      "Get necessary data for interpret call to ensure valid graph for solid state systems and/or "
                      "periodic systems. All necessary data is constructed if not already present. The bond orders are "
                      "constructed with the BondDetector.");
  periodic_system.def("get_data_for_molassembler_interpretation",
                      pybind11::overload_cast<const BondOrderCollection&>(&PeriodicSystem::getDataForMolassemblerInterpretation),
                      pybind11::arg("bond_order_collection"),
                      "Get necessary data for interpret call to ensure valid graph for solid state systems and/or "
                      "periodic systems. Bonds across periodic boundaries have to be negative.");

  periodic_system.def(
      "construct_bond_orders", &PeriodicSystem::constructBondOrders, pybind11::arg("periodic") = true,
      "Constructs bond orders of the atoms without images with negative bond orders across periodic boundaries based "
      "on the BondDetector. The boolean controls whether periodic boundary conditions should be considered for the "
      "bond orders. Bonds across periodic boundaries receive a negative bond order.");

  periodic_system.def(
      "make_bond_orders_across_boundaries_negative", &PeriodicSystem::makeBondOrdersAcrossBoundariesNegative,
      pybind11::arg("bond_order_collection"),
      "Takes bond orders, turns them to absolute values and sets all bond orders to negative values, if the bond "
      "is spanning across periodic boundaries in the current state of the PeriodicSystem");

  // utility functions
  periodic_system.def("center_and_translate_atoms_into_cell", &PeriodicSystem::centerAndTranslateAtomsIntoCell,
                      "First translates center of mass into the center of the cell and then projects all overhanging "
                      "atoms back into the cell");

  // std::vector-like functions
  periodic_system.def("clear", &PeriodicSystem::clear, "Remove all atoms from the collection");

  // Operators
  periodic_system.def(pybind11::self == pybind11::self);
  periodic_system.def(pybind11::self != pybind11::self);
  periodic_system.def(
      "__mul__", [&](const PeriodicSystem& ps, int x) { return ps * x; }, pybind11::is_operator());
  periodic_system.def(
      "__imul__", [&](PeriodicSystem& ps, int x) { return ps *= x; }, pybind11::is_operator());
  periodic_system.def(
      "__mul__", [&](const PeriodicSystem& ps, const Eigen::Vector3i& x) { return ps * x; }, pybind11::is_operator());
  periodic_system.def(
      "__imul__", [&](PeriodicSystem& ps, const Eigen::Vector3i& x) { return ps *= x; }, pybind11::is_operator());

  // Copy support
  periodic_system.def("__copy__", [](const PeriodicSystem& c) -> PeriodicSystem { return PeriodicSystem(c); });
  periodic_system.def("__deepcopy__", [](const PeriodicSystem& c, pybind11::dict /* memo */) -> PeriodicSystem {
    return PeriodicSystem(c);
  });
}
