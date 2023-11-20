/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/DataStructures/SingleParticleEnergies.h"
#include <Utils/Pybind.h>

void init_single_particle_energies(pybind11::module& m) {
  using namespace Scine::Utils;
  pybind11::class_<SingleParticleEnergies> single_particle_energies(m, "SingleParticleEnergies");
  single_particle_energies.def_static("make_unrestricted", &SingleParticleEnergies::createEmptyUnrestrictedEnergies);
  single_particle_energies.def_static("make_restricted", &SingleParticleEnergies::createEmptyRestrictedEnergies);
  single_particle_energies.def_property_readonly("is_restricted", &SingleParticleEnergies::isRestricted);
  single_particle_energies.def_property_readonly("restricted_levels", &SingleParticleEnergies::getRestrictedNLevels);
  single_particle_energies.def_property_readonly("unrestricted_levels", &SingleParticleEnergies::getUnrestrictedNLevels);
  single_particle_energies.def("set_restricted", &SingleParticleEnergies::setRestricted, pybind11::arg("values"));
  single_particle_energies.def("set_unrestricted", &SingleParticleEnergies::setUnrestricted, pybind11::arg("alpha"),
                               pybind11::arg("beta"));
  single_particle_energies.def_property_readonly("restricted_energies", &SingleParticleEnergies::getRestrictedEnergies);
  single_particle_energies.def_property_readonly("alpha", &SingleParticleEnergies::getAlphaEnergies);
  single_particle_energies.def_property_readonly("beta", &SingleParticleEnergies::getBetaEnergies);
}
