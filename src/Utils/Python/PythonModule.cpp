/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Includes */
#include <pybind11/pybind11.h>

void init_adiabatic_mode_localizer(pybind11::module& m);
void init_ao_to_atom_mapping(pybind11::module& m);
void init_atom(pybind11::module& m);
void init_atom_collection(pybind11::module& m);
void init_atomic_gtos(pybind11::module& m);
void init_bond_detector(pybind11::module& m);
void init_bond_order_collection(pybind11::module& m);
void init_conceptual_dft(pybind11::module& m);
void init_constants(pybind11::module& m);
void init_core_submodule(pybind11::module& m);
void init_density_matrix(pybind11::module& m);
void init_element_info(pybind11::module& m);
void init_element_type(pybind11::module& m);
void init_formula_generator(pybind11::module& m);
void init_geometry(pybind11::module& m);
void init_geometry_optimize(pybind11::module& m);
void init_gtf(pybind11::module& m);
void init_gto_expansion(pybind11::module& m);
void init_io(pybind11::module& m);
void init_molecular_orbitals(pybind11::module& m);
void init_molecular_trajectory(pybind11::module& m);
void init_normal_modes(pybind11::module& m);
void init_property_list(pybind11::module& m);
void init_quaternion_fit(pybind11::module& m);
void init_results(pybind11::module& m);
void init_settings(pybind11::module& m);
void init_solute_solvent_complex(pybind11::module& m);
void init_spin_adapted_electronic_transition_result(pybind11::module& m);
void init_spin_adapted_matrix(pybind11::module& m);
void init_structural_completion(pybind11::module& m);
void init_typenames(pybind11::module& m);
void init_thermochemical_calculator(pybind11::module& m);
void init_value_collection(pybind11::module& m);

PYBIND11_MODULE(scine_utilities, m) {
  m.doc() = "Pybind11 Bindings for SCINE-Utils";

  // Cannot follow alphabetical ordering entirely, element types must come first
  init_element_type(m);

  init_adiabatic_mode_localizer(m);
  init_ao_to_atom_mapping(m);
  init_atom(m);
  init_atom_collection(m);
  init_atomic_gtos(m);
  init_bond_detector(m);
  init_bond_order_collection(m);
  init_conceptual_dft(m);
  init_constants(m);
  init_core_submodule(m);
  init_density_matrix(m);
  init_element_info(m);
  init_formula_generator(m);
  init_geometry(m);
  init_gtf(m);
  init_gto_expansion(m);
  init_molecular_orbitals(m);
  init_molecular_trajectory(m);
  init_normal_modes(m);
  init_property_list(m);
  init_quaternion_fit(m);
  init_results(m);
  init_solute_solvent_complex(m);
  init_spin_adapted_electronic_transition_result(m);
  init_spin_adapted_matrix(m);
  init_structural_completion(m);
  init_thermochemical_calculator(m);
  init_typenames(m);
  init_value_collection(m);

  // dependent on several other classes
  init_settings(m);
  init_geometry_optimize(m);
  init_io(m);
}
