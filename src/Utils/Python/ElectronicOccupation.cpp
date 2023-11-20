/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Scf/LcaoUtils/ElectronicOccupation.h"
#include <Utils/Pybind.h>

void init_electronic_occupation(pybind11::module& m) {
  using namespace Scine::Utils;
  using namespace LcaoUtils;

  pybind11::class_<ElectronicOccupation> electronic_occupation(m, "ElectronicOccupation");
  electronic_occupation.def_property_readonly("n_restricted_electrons", &ElectronicOccupation::numberRestrictedElectrons);
  electronic_occupation.def_property_readonly("n_occupied_restricted_orbitals",
                                              &ElectronicOccupation::numberOccupiedRestrictedOrbitals);
  electronic_occupation.def_property_readonly("n_alpha", &ElectronicOccupation::numberAlphaElectrons);
  electronic_occupation.def_property_readonly("n_beta", &ElectronicOccupation::numberBetaElectrons);
  electronic_occupation.def_property_readonly("filled_from_bottom", &ElectronicOccupation::isFilledUpFromTheBottom);
  electronic_occupation.def_property_readonly("has_unpaired_rhf_electron", &ElectronicOccupation::hasUnpairedRHFElectron);
  electronic_occupation.def_property_readonly("has_unpaired_rhf_electron", &ElectronicOccupation::hasUnpairedRHFElectron);
  electronic_occupation.def_property_readonly("restricted", &ElectronicOccupation::isRestricted);
  electronic_occupation.def_property_readonly("unrestricted", &ElectronicOccupation::isUnrestricted);
  electronic_occupation.def("make_unrestricted", &ElectronicOccupation::makeUnrestricted,
                            "From a restricted occupation, transform to an unrestricted occupation");
  electronic_occupation.def("to_unrestricted", &ElectronicOccupation::toUnrestricted,
                            "From a restricted occupation, generate a new unrestricted occupation");
  electronic_occupation.def_property_readonly("filled_restricted_orbitals", &ElectronicOccupation::getFilledRestrictedOrbitals);
  electronic_occupation.def_property_readonly("filled_alpha_orbitals", &ElectronicOccupation::getFilledAlphaOrbitals);
  electronic_occupation.def_property_readonly("filled_beta_orbitals", &ElectronicOccupation::getFilledBetaOrbitals);
  electronic_occupation.def("fill_lowest_unrestricted", &ElectronicOccupation::fillLowestRestrictedOrbitalsWithElectrons,
                            pybind11::arg("n_electrons"));
  electronic_occupation.def("fill_lowest_restricted", &ElectronicOccupation::fillLowestUnrestrictedOrbitals,
                            pybind11::arg("n_alpha"), pybind11::arg("n_beta"));
  electronic_occupation.def("fill_restricted", &ElectronicOccupation::fillSpecifiedRestrictedOrbitals,
                            pybind11::arg("orbitals"));
  electronic_occupation.def("fill_unrestricted", &ElectronicOccupation::fillSpecifiedUnrestrictedOrbitals,
                            pybind11::arg("alpha_orbitals"), pybind11::arg("beta_orbitals"));
}
