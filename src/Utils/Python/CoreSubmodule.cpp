/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <pybind11/pybind11.h>

void init_calculator(pybind11::module& m);
void init_module_manager(pybind11::module& m);
void init_log(pybind11::module& m);

void init_core_submodule(pybind11::module& m) {
  pybind11::module core = m.def_submodule("core");
  core.doc() = R"(
    The ``core`` submodule defines interfaces that permeate the SCINE Project.
    Individual components may offer models of these interfaces that are then
    available through the :class:`core.ModuleManager` class.
  )";

  init_calculator(core);
  init_module_manager(core);
  init_log(core);
}
