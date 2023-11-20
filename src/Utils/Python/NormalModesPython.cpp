/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/GeometricDerivatives/NormalMode.h>
#include <Utils/GeometricDerivatives/NormalModeAnalysis.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_normal_modes(pybind11::module& m) {
  auto normal_modes_submodule = m.def_submodule("normal_modes");

  pybind11::class_<NormalModesContainer> normalModesContainer(normal_modes_submodule, "container");

  normalModesContainer.def("size", &NormalModesContainer::size);
  normalModesContainer.def("get_mode", &NormalModesContainer::getMode, pybind11::arg("mode_index"),
                           "Returns the vibrational mode with index mode_index.");
  normalModesContainer.def("get_mode_as_molecular_trajectory", &NormalModesContainer::getModeAsMolecularTrajectory,
                           pybind11::arg("mode_index"), pybind11::arg("structure"), pybind11::arg("scaling_factor"),
                           "The molecular trajectory representing the mode.");

  normalModesContainer.def("get_wave_numbers", &NormalModesContainer::getWaveNumbers,
                           "Get the wave numbers corresponding to the vibrational modes [cm^(-1)].");

  normal_modes_submodule.def(
      "calculate",
      pybind11::overload_cast<const HessianMatrix&, const AtomCollection&>(&NormalModeAnalysis::calculateNormalModes),
      pybind11::arg("hessian"), pybind11::arg("atoms"), "Calculate the mass weighted normal modes.");
  normal_modes_submodule.def(
      "calculate",
      pybind11::overload_cast<const HessianMatrix&, const ElementTypeCollection&, const PositionCollection&, bool>(
          &NormalModeAnalysis::calculateNormalModes),
      pybind11::arg("hessian"), pybind11::arg("elements"), pybind11::arg("positions"),
      pybind11::arg("normalize") = false, "Calculate the mass weighted normal modes.");
  normal_modes_submodule.def(
      "calculate",
      pybind11::overload_cast<const PartialHessian&, const AtomCollection&>(&NormalModeAnalysis::calculateNormalModes),
      pybind11::arg("partial_hessian"), pybind11::arg("atoms"),
      "Calculate the mass weighted normal modes from a"
      "partial Hessian and the super system.");
  normal_modes_submodule.def(
      "calculate",
      pybind11::overload_cast<const PartialHessian&, const ElementTypeCollection&, const PositionCollection&, bool>(
          &NormalModeAnalysis::calculateNormalModes),
      pybind11::arg("partial_hessian"), pybind11::arg("elements"), pybind11::arg("positions"),
      pybind11::arg("normalize") = false,
      "Calculate the mass weighted normal modes from a partial Hessian and the super system.");

  normal_modes_submodule.def("get_harmonic_inversion_point", &NormalModeAnalysis::calculateHarmonicInversionPoint,
                             pybind11::arg("wavenumber"), pybind11::arg("n"),
                             "Returns the n-th harmonic inversion point displacement for a given wavenumber.");

  pybind11::class_<NormalMode> normalMode(normal_modes_submodule, "mode");

  normalMode.def(pybind11::init<double, DisplacementCollection>(), pybind11::arg("wave_number"), pybind11::arg("mode"),
                 "Initialize a new normal mode object.");
  normalMode.def("get_wave_number", &NormalMode::getWaveNumber, "Get the wave number.");
  normalMode.def("get_mode", &NormalMode::getMode, "Get mode.");
}
