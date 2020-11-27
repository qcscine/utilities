/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Math/IterativeDiagonalizer/SpinAdaptedEigenContainer.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_spin_adapted_electronic_transition_result(pybind11::module& m) {
  pybind11::class_<ElectronicTransitionResult> electronicTransitionResult(m, "ElectronicTransitionResult");
  electronicTransitionResult.def_readonly("eigenstates", &ElectronicTransitionResult::eigenStates);
  electronicTransitionResult.def_readonly("transition_dipoles", &ElectronicTransitionResult::transitionDipoles);

  pybind11::class_<SpinAdaptedElectronicTransitionResult> spinAdaptedElectronicTransitionResult(
      m, "SpinAdaptedElectronicTransitionResult");
  spinAdaptedElectronicTransitionResult.def_readonly("singlet", &SpinAdaptedElectronicTransitionResult::singlet);
  spinAdaptedElectronicTransitionResult.def_readonly("triplet", &SpinAdaptedElectronicTransitionResult::triplet);
}