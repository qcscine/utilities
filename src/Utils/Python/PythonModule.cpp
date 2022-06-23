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
void init_atomic_second_derivative_collection(pybind11::module& m);
void init_bond_detector(pybind11::module& m);
void init_bond_order_collection(pybind11::module& m);
void init_bspline_functionalities(pybind11::module& m);
void init_conceptual_dft(pybind11::module& m);
void init_constants(pybind11::module& m);
void init_core_submodule(pybind11::module& m);
void init_cp2k_external_qc(pybind11::module& m);
void init_density_matrix(pybind11::module& m);
void init_descriptor_collection(pybind11::module& m);
void init_dipole_matrix(pybind11::module& m);
void init_electronic_occupation(pybind11::module& m);
void init_element_info(pybind11::module& m);
void init_element_type(pybind11::module& m);
void init_formula_generator(pybind11::module& m);
void init_geometry(pybind11::module& m);
void init_geometry_optimize(pybind11::module& m);
void init_gtf(pybind11::module& m);
void init_gto_expansion(pybind11::module& m);
void init_io(pybind11::module& m);
void init_iterative_diagonalizer(pybind11::module& m);
void init_molecular_dynamics(pybind11::module& m);
void init_molecular_orbitals(pybind11::module& m);
void init_molecular_surface(pybind11::module& m);
void init_molecular_trajectory(pybind11::module& m);
void init_normal_modes(pybind11::module& m);
void init_orca_parser(pybind11::module& m);
void init_periodic_boundaries(pybind11::module& m);
void init_periodic_system(pybind11::module& m);
void init_property_list(pybind11::module& m);
void init_quaternion_fit(pybind11::module& m);
void init_results(pybind11::module& m);
void init_setting_descriptors(pybind11::module& m);
void init_settings(pybind11::module& m);
void init_settings_names(pybind11::module& m);
void init_single_particle_energies(pybind11::module& m);
void init_solute_solvent_complex(pybind11::module& m);
void init_spin_adapted_electronic_transition_result(pybind11::module& m);
void init_spin_adapted_matrix(pybind11::module& m);
void init_structural_completion(pybind11::module& m);
void init_thermochemical_calculator(pybind11::module& m);
void init_typenames(pybind11::module& m);
void init_value_collection(pybind11::module& m);
void init_regression_functionalities(pybind11::module& m);

PYBIND11_MODULE(scine_utilities, m) {
  m.doc() = "Pybind11 Bindings for SCINE-Utilities";

  // Ordering is important!
  init_element_type(m);
  init_ao_to_atom_mapping(m);
  init_atom(m);
  init_atom_collection(m);
  init_bond_order_collection(m);
  init_conceptual_dft(m);
  init_constants(m);
  init_density_matrix(m);
  init_dipole_matrix(m);
  init_element_info(m);
  init_formula_generator(m);
  init_molecular_trajectory(m);
  init_periodic_boundaries(m);
  init_periodic_system(m);
  init_geometry(m);
  init_gtf(m);
  init_gto_expansion(m);
  init_atomic_gtos(m);
  init_molecular_orbitals(m);
  init_normal_modes(m);
  init_bond_detector(m);
  init_property_list(m);
  init_quaternion_fit(m);
  init_typenames(m);
  init_thermochemical_calculator(m);
  init_spin_adapted_electronic_transition_result(m);
  init_spin_adapted_matrix(m);
  init_single_particle_energies(m);
  init_orca_parser(m);
  init_electronic_occupation(m);
  init_atomic_second_derivative_collection(m);
  init_results(m);
  init_structural_completion(m);
  init_value_collection(m);
  init_setting_descriptors(m);
  init_settings(m);
  init_settings_names(m);
  init_core_submodule(m);
  init_cp2k_external_qc(m);
  init_molecular_dynamics(m);
  init_molecular_surface(m);
  init_solute_solvent_complex(m);
  init_adiabatic_mode_localizer(m);
  init_geometry_optimize(m);
  init_io(m);
  init_iterative_diagonalizer(m);
  init_bspline_functionalities(m);
  init_regression_functionalities(m);
}
