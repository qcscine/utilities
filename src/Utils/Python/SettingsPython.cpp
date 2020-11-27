/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Settings.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

extern void update(UniversalSettings::ValueCollection& coll, const pybind11::dict& dictionary, bool preserveTypes);

namespace {

Settings fromDict(const std::string& name, const pybind11::dict& dict) {
  Settings settings(name);
  update(settings, dict, false);
  return settings;
}

} // namespace

void init_settings(pybind11::module& m) {
  pybind11::class_<Settings, UniversalSettings::ValueCollection, std::shared_ptr<Settings>> settings(m, "Settings",
                                                                                                     R"delim(
      A named :class:`ValueCollection`, with configurable key name and value or
      type restrictions on entries.

      >>> settings = Settings("example_settings")
      >>> settings.name
      'example_settings'
      >>> settings["molecular_charge"] = 4
      >>> settings["spin_multiplicity"] = 1
      >>> settings.update({"molecular_charge": 0})
      >>> settings["molecular_charge"]
      0
    )delim");

  settings.def(pybind11::init<std::string>(), pybind11::arg("name"), "Initialize a named settings object");

  settings.def(pybind11::init(&fromDict), pybind11::arg("name"), pybind11::arg("dict"),
               "Initialize settings from name and a dictionary of values");

  settings.def_property_readonly("name", &Settings::name, "Name of the settings");

  settings.def("check", &Settings::check, "Checks if the settings are acceptable w.r.t. the defined boundaries");
}
