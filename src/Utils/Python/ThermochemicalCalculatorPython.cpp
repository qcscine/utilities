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
  pybind11::class_<ThermochemicalContainer> thermochemicalContainer(m, "ThermochemicalContainer");

  thermochemicalContainer.def_readonly("entropy", &ThermochemicalContainer::entropy);
  thermochemicalContainer.def_readonly("enthalpy", &ThermochemicalContainer::enthalpy);
  thermochemicalContainer.def_readonly("heat_capacity_p", &ThermochemicalContainer::heatCapacityP);
  thermochemicalContainer.def_readonly("heat_capacity_v", &ThermochemicalContainer::heatCapacityV);
  thermochemicalContainer.def_readonly("gibbs_free_energy", &ThermochemicalContainer::gibbsFreeEnergy);
  thermochemicalContainer.def_readonly("zero_point_vibrational_energy", &ThermochemicalContainer::zeroPointVibrationalEnergy);

  pybind11::class_<ThermochemicalComponentsContainer> thermochemicalComponentsContainer(
      m, "ThermochemicalComponentsContainer");

  thermochemicalComponentsContainer.def_readonly("vibrational_component",
                                                 &ThermochemicalComponentsContainer::vibrationalComponent);
  thermochemicalComponentsContainer.def_readonly("rotational_component", &ThermochemicalComponentsContainer::rotationalComponent);
  thermochemicalComponentsContainer.def_readonly("translational_component",
                                                 &ThermochemicalComponentsContainer::translationalComponent);
  thermochemicalComponentsContainer.def_readonly("electronic_component", &ThermochemicalComponentsContainer::electronicComponent);
  thermochemicalComponentsContainer.def_readonly("overall", &ThermochemicalComponentsContainer::overall);
}
