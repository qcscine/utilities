/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/DataStructures/PeriodicBoundaries.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Geometry/PeriodicSystem.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_periodic_system(pybind11::module& m) {
  pybind11::class_<PeriodicSystem> periodic_system(m, "PeriodicSystem",
                                                   R"delim(
      A class representing a collection of Atoms including periodic boundary conditions.
      Holds the AtomCollection, PeriodicBoundaries, and the set of indices representing solid state atom indices
      as public members.
      Additionally, includes functionalities based on periodic boundaries focused on structures such as
      primitive cell reduction, generation of image atoms for graph consistency (also a call that allows to directly
      generate the necessary data for SCINE Molassembler), and comparison and supercell operations.

      All method calls that may generate a new bond order collection also allow to detect bonds between
      solid state atoms based on van der Waals radii instead of nearest neighbors (optional flag, true per default).
      This will likely overestimate the bonding within a solid state structure, but also avoids that solid state
      structures may be split up into separate graphs.
      The van der Waals radii are specified in scine_utilities.ElementInfo
    )delim");
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
                      pybind11::arg("use_solid_state_van_der_waals_bonds") = true, "Get the atoms plus the image atoms");
  periodic_system.def("get_image_atoms", &PeriodicSystem::getImageAtoms,
                      pybind11::arg("use_solid_state_van_der_waals_bonds") = true, "Get the image atoms only");
  periodic_system.def("get_image_atoms_map", &PeriodicSystem::getImageAtomsMap,
                      pybind11::arg("use_solid_state_van_der_waals_bonds") = true,
                      "Get a map pointing from the index of the image atom to the index of the real space atom");

  // methods surrounding image atoms
  periodic_system.def("get_data_for_molassembler_interpretation",
                      pybind11::overload_cast<bool>(&PeriodicSystem::getDataForMolassemblerInterpretation),
                      pybind11::arg("use_solid_state_van_der_waals_bonds") = true,
                      "Get necessary data for interpret call to ensure valid graph for solid state systems and/or "
                      "periodic systems. All necessary data is constructed if not already present. The bond orders are "
                      "constructed with the BondDetector.");
  periodic_system.def(
      "get_data_for_molassembler_interpretation",
      pybind11::overload_cast<const BondOrderCollection&, bool>(&PeriodicSystem::getDataForMolassemblerInterpretation),
      pybind11::arg("bond_order_collection"), pybind11::arg("remove_solid_second_shell") = false,
      "Get necessary data for interpret call to ensure valid graph for solid state systems and/or "
      "periodic systems. Bonds across periodic boundaries have to be negative.");

  periodic_system.def(
      "construct_bond_orders", &PeriodicSystem::constructBondOrders, pybind11::arg("periodic") = true,
      pybind11::arg("use_solid_state_van_der_waals_bonds") = true,
      "Constructs bond orders of the atoms without images with negative bond orders across periodic boundaries based "
      "on the BondDetector. The first boolean controls whether periodic boundary conditions should be considered for "
      "the "
      "bond orders. The second boolean controls if bonds between solid state atoms should be evaluated by"
      "nearest neighbors or van der Waals radii. Bonds across periodic boundaries receive a negative bond order.");

  periodic_system.def(
      "make_bond_orders_across_boundaries_negative", &PeriodicSystem::makeBondOrdersAcrossBoundariesNegative,
      pybind11::arg("bond_order_collection"),
      "Takes bond orders, turns them to absolute values and sets all bond orders to negative values, if the bond "
      "is spanning across periodic boundaries in the current state of the PeriodicSystem");

  periodic_system.def("get_primitive_cell_system", &PeriodicSystem::getPrimitiveCellSystem, pybind11::arg("epsilon") = 1e-6,
                      pybind11::arg("solid_state_only") = false, "Get the system reduced to the primitive cell.");

  // utility functions
  periodic_system.def("center_and_translate_atoms_into_cell", &PeriodicSystem::centerAndTranslateAtomsIntoCell,
                      "First translates center of mass into the center of the cell and then projects all overhanging "
                      "atoms back into the cell");

  // std::vector-like functions
  periodic_system.def("clear", &PeriodicSystem::clear, "Remove all atoms from the collection");

  // Operators
  periodic_system.def("is_approx", &PeriodicSystem::isApprox, pybind11::arg("other_system"),
                      pybind11::arg("epsilon") = 1e-6, "Allows to set the accuracy of the fuzzy comparison.");
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
  periodic_system.def(
      "__mul__", [&](const PeriodicSystem& ps, const std::vector<int>& x) { return ps * x; }, pybind11::is_operator());
  periodic_system.def(
      "__imul__", [&](PeriodicSystem& ps, const std::vector<int>& x) { return ps *= x; }, pybind11::is_operator());

  // Copy support
  periodic_system.def("__copy__", [](const PeriodicSystem& c) -> PeriodicSystem { return PeriodicSystem(c); });
  periodic_system.def("__deepcopy__", [](const PeriodicSystem& c, pybind11::dict /* memo */) -> PeriodicSystem {
    return PeriodicSystem(c);
  });
  periodic_system.def(pybind11::pickle(
      [](const PeriodicSystem& ps) { // __getstate__
        /* Return a tuple that fully encodes the state of the object */
        std::vector<std::string> elements;
        for (const auto& e : ps.atoms.getElements()) {
          elements.push_back(ElementInfo::symbol(e));
        }
        return pybind11::make_tuple(ps.pbc.getCellMatrix(), elements, ps.atoms.getPositions(), ps.solidStateAtomIndices);
      },
      [](pybind11::tuple t) { // __setstate__
        if (t.size() != 4)
          throw std::runtime_error("Invalid state for PeriodicSystem!");
        ElementTypeCollection elements;
        for (const auto& e : t[1].cast<std::vector<std::string>>()) {
          elements.push_back(ElementInfo::elementTypeForSymbol(e));
        }
        PeriodicSystem ps(PeriodicBoundaries(t[0].cast<Eigen::Matrix3d>()), elements, t[2].cast<PositionCollection>(),
                          t[3].cast<std::unordered_set<unsigned>>());
        return ps;
      }));
}
