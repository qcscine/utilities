/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Includes */
#include <pybind11/pybind11.h>

void init_atom(pybind11::module& m);
void init_atom_collection(pybind11::module& m);
void init_bond_detector(pybind11::module& m);
void init_bond_order_collection(pybind11::module& m);
void init_constants(pybind11::module& m);
void init_element_data(pybind11::module& m);
void init_element_info(pybind11::module& m);
void init_element_type(pybind11::module& m);
void init_formula_generator(pybind11::module& m);
void init_geometry(pybind11::module& m);
void init_io(pybind11::module& m);
void init_logger(pybind11::module& m);
void init_molecular_trajectory(pybind11::module& m);
void init_property_list(pybind11::module& m);
void init_results(pybind11::module& m);
void init_structural_completion(pybind11::module& m);
void init_typenames(pybind11::module& m);
void init_core_submodule(pybind11::module& m);
void init_thermochemical_calculator(pybind11::module& m);
void init_spin_adapted_matrix(pybind11::module& m);
void init_spin_adapted_electronic_transition_result(pybind11::module& m);

PYBIND11_MODULE(scine_utilities, m) {
  m.doc() = "Pybind11 Bindings for SCINE-Utils";

  // Cannot follow alphabetical ordering entirely, element types must come first
  init_element_type(m);

  init_atom(m);
  init_atom_collection(m);
  init_bond_detector(m);
  init_bond_order_collection(m);
  init_constants(m);
  init_core_submodule(m);
  init_element_data(m);
  init_element_info(m);
  init_formula_generator(m);
  init_geometry(m);
  init_logger(m);
  init_molecular_trajectory(m);
  init_spin_adapted_matrix(m);
  init_thermochemical_calculator(m);
  init_spin_adapted_electronic_transition_result(m);
  init_property_list(m);
  init_results(m);
  init_structural_completion(m);
  init_typenames(m);

  // io is dependent on several other classes
  init_io(m);
}
