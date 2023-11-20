/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/AtomicGtos.h>
#include <Utils/DataStructures/Gtf.h>
#include <Utils/DataStructures/GtoExpansion.h>
#include <Utils/Pybind.h>
#include <pybind11/complex.h>

using namespace Scine::Utils;

void init_atomic_gtos(pybind11::module& m) {
  pybind11::class_<AtomicGtos> atomic_gtos(m, "AtomicGtos");

  atomic_gtos.def_readonly(
      "s", &AtomicGtos::s,
      "Optional type allows to check for existence of the individual expansion with if. However, there is a bug in "
      "pybind11 version <= 2.5, which causes the gtfs to be deleted when accessing any expansion member");

  atomic_gtos.def_readonly(
      "p", &AtomicGtos::p,
      "Optional type allows to check for existence of the individual expansion with if. However, there is a bug in "
      "pybind11 version <= 2.5, which causes the gtfs to be deleted when accessing any expansion member");
  atomic_gtos.def_readonly(
      "d", &AtomicGtos::d,
      "Optional type allows to check for existence of the individual expansion with if. However, there is a bug in "
      "pybind11 version <= 2.5, which causes the gtfs to be deleted when accessing any expansion member");

  atomic_gtos.def("get_nwchem_format", &AtomicGtos::getNWChemFormat,
                  "get basis in fitting NWChem format for atomic_gto of one element, but still needs to be combined "
                  "into dictionary with elements as keys.");
  atomic_gtos.def("get_gtfs", &AtomicGtos::getGtfs, "Returns dictionary of expansions with list of Gtfs.");
}

void init_gto_expansion(pybind11::module& m) {
  pybind11::class_<GtoExpansion> gto_expansion(m, "GtoExpansion");

  gto_expansion.def_readonly("angular_momentum", &GtoExpansion::angularMomentum, "angular momentum as an integer.");
  gto_expansion.def_readonly("gtfs", &GtoExpansion::gtfs,
                             "List of all Gtfs of expansion. Be aware that you need Pybind11 version > 2.5 to be able "
                             "to use this member variable.");
  gto_expansion.def_property_readonly("n_aos", &GtoExpansion::nAOs, "number of spherical AOs");
}

void init_gtf(pybind11::module& m) {
  pybind11::class_<Gtf> gtf(m, "Gtf");

  gtf.def_readonly("exponent", &Gtf::exponent);
  gtf.def_readonly("coefficient", &Gtf::coefficient);
  gtf.def_readonly("normalized_coefficient", &Gtf::normalizedCoefficient);
}
