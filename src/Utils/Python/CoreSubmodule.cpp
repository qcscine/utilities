/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <pybind11/pybind11.h>

void init_log(pybind11::module& m);
void init_object_with_log(pybind11::module& m);
void init_calculator(pybind11::module& m);
void init_calculator_with_reference(pybind11::module& m);
void init_module_manager(pybind11::module& m);

void init_core_submodule(pybind11::module& m) {
  pybind11::module core = m.def_submodule("core");
  core.doc() = R"(
    The ``core`` submodule defines interfaces that permeate the SCINE Project.
    Individual components may offer models of these interfaces that are then
    available through the :class:`core.ModuleManager` class.
  )";

  init_log(core);
  init_calculator(core);
  init_calculator_with_reference(core);
  init_module_manager(core);
}
