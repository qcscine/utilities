/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <pybind11/pybind11.h>

void init_calculator(pybind11::module& m);
void init_module_manager(pybind11::module& m);

void init_core_submodule(pybind11::module& m) {
  pybind11::module core = m.def_submodule("core");

  init_calculator(core);
  init_module_manager(core);
}
