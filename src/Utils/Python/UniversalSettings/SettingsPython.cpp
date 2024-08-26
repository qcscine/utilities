/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Pybind.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/GenericValueVariant.h>

using namespace Scine::Utils;

// Defined in ValueCollectionPython.cpp
extern void update(UniversalSettings::ValueCollection& coll, const pybind11::dict& dictionary, bool preserveTypes);

namespace {

Settings fromDict(const std::string& name, const pybind11::dict& dict) {
  Settings settings(name);
  update(settings, dict, false);
  return settings;
}

} // namespace

void init_settings(pybind11::module& m) {
  using namespace UniversalSettings;
  pybind11::class_<Settings, ValueCollection, std::shared_ptr<Settings>> settings(m, "Settings",
                                                                                  R"delim(
      A :class:`ValueCollection` with restrictions on keys and value types.

      >>> descriptors = DescriptorCollection("example_settings")
      >>> descriptors["molecular_charge"] = IntDescriptor("Charge on a molecule")
      >>> descriptors["spin_multiplicity"] = IntDescriptor("Spin multiplicity")
      >>> descriptors["spin_multiplicity"].minimum = 1
      >>> settings = Settings(descriptors)
      >>> settings["spin_multiplicity"]
      1
      >>> settings["molecular_charge"] = -4
      >>> settings["spin_multiplicity"] = 3
      >>> settings.update({"molecular_charge": 2})  # Only existing fields are updated
      >>> settings["molecular_charge"]
      2
    )delim");

  settings.def(pybind11::init<ValueCollection, DescriptorCollection>(), pybind11::arg("values"), pybind11::arg("fields"));

  settings.def(pybind11::init([](DescriptorCollection descriptors) {
                 auto coll = createDefaultValueCollection(descriptors);
                 return Settings{std::move(coll), std::move(descriptors)};
               }),
               pybind11::arg("fields"));

  settings.def(pybind11::init(&fromDict), pybind11::arg("name"), pybind11::arg("dict"),
               "Initialize settings from name and a dictionary of values");

  settings.def(
      "update",
      [](Settings& settings, const pybind11::dict& dict) {
        ValueCollection coll;
        update(coll, dict, false);
        settings.merge(coll, true);
      },
      R"delim(
      Update existing keys with new values from a dictionary

      ..note: Values without matching key are ignored.
    )delim");

  settings.def(
      "update", [](Settings& settings, const ValueCollection& coll) { settings.merge(coll, true); },
      pybind11::arg("collection"),
      R"delim(
      Update existing keys with new values from a ValueCollection

      ..note: Values without matching key in settings are ignored.
    )delim");

  settings.def_property_readonly("descriptor_collection", &Settings::getDescriptorCollection);

  settings.def("check", &Settings::valid, "Checks if the settings are acceptable w.r.t. the defined boundaries");

  settings.def("valid", &Settings::valid, "Checks if the settings are acceptable w.r.t. the defined boundaries");

  settings.def("throw_incorrect_settings", &Settings::throwIncorrectSettings,
               "Throws an exception with the incorrect setting; settings must be invalid");

  settings.def("reset", &Settings::resetToDefaults, "Resets the settings to the defaults");

  settings.def(
      "extract",
      [](Settings& settings, const std::string& name, const pybind11::object& default_value) -> pybind11::object {
        if (default_value.is(pybind11::none())) {
          if (settings.valueExists(name)) {
            auto result = settings.getValue(name);
            settings.dropValue(name);
            return pybind11::cast(GenericValueMeta::convert(result));
          }
          return default_value;
        }
        GenericValue result =
            settings.extract(name, GenericValueMeta::convert(default_value.cast<GenericValueMeta::Variant>()));
        return pybind11::cast(GenericValueMeta::convert(result));
      },
      pybind11::arg("name"), pybind11::arg("default_value") = pybind11::none(),
      "Gets a value from the settings based on the given key and removes it from the settings.\n"
      "If it is not present the given default value is returned.");

  settings.def("__setitem__", [](Settings& settings, const std::string& key, const GenericValueMeta::Variant& v) {
    if (settings.valueExists(key)) {
      settings.modifyValue(key, GenericValueMeta::convert(v));
    }
    else {
      const auto& descriptors = settings.getDescriptorCollection();
      const auto findIter = std::find_if(std::begin(descriptors), std::end(descriptors),
                                         [&](const auto& keyValuePair) { return keyValuePair.first == key; });

      if (findIter == std::end(descriptors)) {
        throw std::runtime_error("There is no matching key '" + key + "' in these Settings!");
      }

      settings.addGenericValue(key, GenericValueMeta::convert(v));
    }
  });

  settings.def("__repr__", [](pybind11::handle h) -> std::string {
    const std::string name = qualifiedName(h);
    const Settings& settings = h.cast<Settings>();
    if (settings.empty()) {
      std::string repr = name + "('";
      repr += settings.name() + "', {})";
      return repr;
    }

    return "<" + name + " instance>";
  });

  /* Pickling support would be better, but is very tricky seeing as
   * ValueCollections have nested class states, and using C++ copy constructors
   * might be preferable
   */
  settings.def("__deepcopy__",
               [](const Settings& settings, const pybind11::dict & /* memo */) -> Settings { return settings; });
}
