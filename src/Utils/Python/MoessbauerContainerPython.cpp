/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/ExternalQC/Orca/MoessbauerParameters.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/optional.hpp>

void init_moessbauer_container(pybind11::module& m) {
  using namespace Scine::Utils::ExternalQC::Moessbauer;
  pybind11::class_<MoessbauerParameterContainer> moessbauerParameterContainer(m, "MoessbauerParameterContainer");
  moessbauerParameterContainer.def(pybind11::init<>());
  moessbauerParameterContainer.def_readonly("num_irons", &MoessbauerParameterContainer::numIrons);
  moessbauerParameterContainer.def_readonly("quadrupole_splittings", &MoessbauerParameterContainer::quadrupoleSplittings);
  moessbauerParameterContainer.def_readonly("etas", &MoessbauerParameterContainer::etas);
  moessbauerParameterContainer.def_readonly("densities", &MoessbauerParameterContainer::densities);
}
