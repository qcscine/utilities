/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Properties/Thermochemistry/ThermochemistryCalculator.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/optional.hpp>

using namespace Scine::Utils;

void init_thermochemical_calculator(pybind11::module& m) {
  pybind11::class_<ThermochemistryCalculator> thermochemistryCalculator(m, "ThermochemistryCalculator");
  pybind11::class_<ThermochemicalContainer> thermochemicalContainer(m, "ThermochemicalContainer");
  pybind11::class_<ThermochemicalComponentsContainer> thermochemicalComponentsContainer(
      m, "ThermochemicalComponentsContainer");

  thermochemistryCalculator.def(pybind11::init<const HessianMatrix&, const AtomCollection&, int, double>(),
                                pybind11::arg("hessian"), pybind11::arg("atoms"), pybind11::arg("multiplicity"),
                                pybind11::arg("energy"), "Initialize a thermochemistry calculator");
  thermochemistryCalculator.def(
      pybind11::init<const HessianMatrix&, ElementTypeCollection, const PositionCollection&, int, double>(),
      pybind11::arg("hessian"), pybind11::arg("elements"), pybind11::arg("positions"), pybind11::arg("multiplicity"),
      pybind11::arg("energy"), "Initialize a thermochemistry calculator");
  thermochemistryCalculator.def("set_temperature", &ThermochemistryCalculator::setTemperature);
  thermochemistryCalculator.def("set_molecular_symmetry", &ThermochemistryCalculator::setMolecularSymmetryNumber);
  thermochemistryCalculator.def("calculate", &ThermochemistryCalculator::calculate);

  thermochemicalContainer.def_readonly("entropy", &ThermochemicalContainer::entropy);
  thermochemicalContainer.def_readonly("enthalpy", &ThermochemicalContainer::enthalpy);
  thermochemicalContainer.def_readonly("heat_capacity_p", &ThermochemicalContainer::heatCapacityP);
  thermochemicalContainer.def_readonly("heat_capacity_v", &ThermochemicalContainer::heatCapacityV);
  thermochemicalContainer.def_readonly("gibbs_free_energy", &ThermochemicalContainer::gibbsFreeEnergy);
  thermochemicalContainer.def_readonly("zero_point_vibrational_energy", &ThermochemicalContainer::zeroPointVibrationalEnergy);
  thermochemicalContainer.def_readonly("symmetry_number", &ThermochemicalContainer::symmetryNumber);

  thermochemicalComponentsContainer.def_readonly("vibrational_component",
                                                 &ThermochemicalComponentsContainer::vibrationalComponent);
  thermochemicalComponentsContainer.def_readonly("rotational_component", &ThermochemicalComponentsContainer::rotationalComponent);
  thermochemicalComponentsContainer.def_readonly("translational_component",
                                                 &ThermochemicalComponentsContainer::translationalComponent);
  thermochemicalComponentsContainer.def_readonly("electronic_component", &ThermochemicalComponentsContainer::electronicComponent);
  thermochemicalComponentsContainer.def_readonly("overall", &ThermochemicalComponentsContainer::overall);
}
