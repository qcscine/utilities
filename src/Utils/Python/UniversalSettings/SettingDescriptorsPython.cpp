/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Pybind.h>
#include <Utils/UniversalSettings/BoolDescriptor.h>
#include <Utils/UniversalSettings/CollectionListDescriptor.h>
#include <Utils/UniversalSettings/DirectoryDescriptor.h>
#include <Utils/UniversalSettings/DoubleDescriptor.h>
#include <Utils/UniversalSettings/DoubleListDescriptor.h>
#include <Utils/UniversalSettings/FileDescriptor.h>
#include <Utils/UniversalSettings/GenericDescriptorVariant.h>
#include <Utils/UniversalSettings/GenericValueVariant.h>
#include <Utils/UniversalSettings/IntDescriptor.h>
#include <Utils/UniversalSettings/IntListDescriptor.h>
#include <Utils/UniversalSettings/OptionListDescriptor.h>
#include <Utils/UniversalSettings/ParametrizedOptionListDescriptor.h>
#include <Utils/UniversalSettings/SettingDescriptor.h>
#include <Utils/UniversalSettings/StringDescriptor.h>
#include <Utils/UniversalSettings/StringListDescriptor.h>

using namespace Scine::Utils;

using namespace UniversalSettings;

namespace {

void init_setting_descriptor(pybind11::module& m) {
  pybind11::class_<SettingDescriptor, std::shared_ptr<SettingDescriptor>> settingDescriptor(
      m, "SettingDescriptor", "Base class to represent a setting's type and possible values");

  settingDescriptor.def_property("property_description", &SettingDescriptor::getPropertyDescription,
                                 &SettingDescriptor::setPropertyDescription, "Explanation of what the setting is for");

  settingDescriptor.def_property_readonly("default_generic_value", [](const SettingDescriptor& descriptor) {
    return GenericValueMeta::convert(descriptor.getDefaultGenericValue());
  });

  settingDescriptor.def(
      "valid_value",
      [](const SettingDescriptor& descriptor, const GenericValueMeta::Variant& variant) {
        return descriptor.validValue(GenericValueMeta::convert(variant));
      },
      pybind11::arg("value"), "Checks if a particular type-erased value is valid for this setting");
}

void init_fundamental_descriptors(pybind11::module& m) {
  pybind11::class_<BoolDescriptor, SettingDescriptor, std::shared_ptr<BoolDescriptor>> boolDescriptor(m,
                                                                                                      "BoolDescriptor");

  boolDescriptor.def(pybind11::init<std::string>(), pybind11::arg("description"));

  boolDescriptor.def_property("default_value", &BoolDescriptor::getDefaultValue, &BoolDescriptor::setDefaultValue,
                              "Default value of the setting");

  pybind11::class_<IntDescriptor, SettingDescriptor, std::shared_ptr<IntDescriptor>> intDescriptor(m, "IntDescriptor");

  intDescriptor.def(pybind11::init<std::string>(), pybind11::arg("description"));

  intDescriptor.def_property("default_value", &IntDescriptor::getDefaultValue, &IntDescriptor::setDefaultValue,
                             "Default value of the setting");

  intDescriptor.def_property("minimum", &IntDescriptor::getMinimum, &IntDescriptor::setMinimum, "Lower bound on valid values");

  intDescriptor.def_property("maximum", &IntDescriptor::getMaximum, &IntDescriptor::setMaximum, "Upper bound on valid values");

  intDescriptor.def("valid_value", pybind11::overload_cast<int>(&IntDescriptor::validValue, pybind11::const_));

  pybind11::class_<DoubleDescriptor, SettingDescriptor, std::shared_ptr<DoubleDescriptor>> doubleDescriptor(
      m, "DoubleDescriptor");

  doubleDescriptor.def(pybind11::init<std::string>(), pybind11::arg("description"));

  doubleDescriptor.def_property("default_value", &DoubleDescriptor::getDefaultValue, &DoubleDescriptor::setDefaultValue,
                                "Default value of the settings");

  doubleDescriptor.def_property("minimum", &DoubleDescriptor::getMinimum, &DoubleDescriptor::setMinimum,
                                "Lower bound on valid values");

  doubleDescriptor.def_property("maximum", &DoubleDescriptor::getMaximum, &DoubleDescriptor::setMaximum,
                                "Upper bound on valid values");

  doubleDescriptor.def("valid_value", pybind11::overload_cast<double>(&DoubleDescriptor::validValue, pybind11::const_));
}

void init_string_descriptor(pybind11::module& m) {
  pybind11::class_<StringDescriptor, SettingDescriptor, std::shared_ptr<StringDescriptor>> stringDescriptor(
      m, "StringDescriptor");

  stringDescriptor.def(pybind11::init<std::string>(), pybind11::arg("description"));

  stringDescriptor.def_property("default_value", &StringDescriptor::getDefaultValue, &StringDescriptor::setDefaultValue,
                                "Default value of the settings");
}

void init_file_descriptor(pybind11::module& m) {
  pybind11::class_<FileDescriptor, SettingDescriptor, std::shared_ptr<FileDescriptor>> fileDescriptor(m,
                                                                                                      "FileDescriptor");

  fileDescriptor.def(pybind11::init<std::string>(), pybind11::arg("description"));

  fileDescriptor.def_property("default_value", &FileDescriptor::getDefaultValue, &FileDescriptor::setDefaultValue,
                              "Default value of the settings");

  fileDescriptor.def_property("file_must_already_exist", &FileDescriptor::fileMustAlreadyExist,
                              &FileDescriptor::setFileMustAlreadyExist, "Whether the referenced file must already exist");

  pybind11::enum_<FileDescriptor::FileType> fileType(fileDescriptor, "FileType");
  fileType.value("Any", FileDescriptor::FileType::Any);
  fileType.value("Executable", FileDescriptor::FileType::Executable);

  fileDescriptor.def_property("file_type", &FileDescriptor::getFileType, &FileDescriptor::setFileType,
                              "What kind of file is referenced");

  fileDescriptor.def_property_readonly("name_filters", &FileDescriptor::getNameFilters,
                                       "Qt5 name filter strings, used for GUI only");

  fileDescriptor.def("add_name_filter", &FileDescriptor::addNameFilter, "Add a Qt5 name filter string");
}

void init_directory_descriptor(pybind11::module& m) {
  pybind11::class_<DirectoryDescriptor, SettingDescriptor, std::shared_ptr<DirectoryDescriptor>> directoryDescriptor(
      m, "DirectoryDescriptor");

  directoryDescriptor.def(pybind11::init<std::string>(), pybind11::arg("description"));
  directoryDescriptor.def_property("default_value", &DirectoryDescriptor::getDefaultValue,
                                   &DirectoryDescriptor::setDefaultValue, "Default value of the settings");
}

void init_option_list_descriptor(pybind11::module& m) {
  pybind11::class_<OptionListDescriptor, SettingDescriptor, std::shared_ptr<OptionListDescriptor>> optionListDescriptor(
      m, "OptionListDescriptor");

  optionListDescriptor.def(pybind11::init<std::string>(), pybind11::arg("description"));

  optionListDescriptor.def_property("default_value", &OptionListDescriptor::getDefaultOption,
                                    &OptionListDescriptor::setDefaultOption, "Default value of the settings");

  optionListDescriptor.def_property_readonly("options", &OptionListDescriptor::getAllOptions);

  optionListDescriptor.def("add_option", &OptionListDescriptor::addOption, pybind11::arg("option"));
}

void init_parametrized_option_list_descriptor(pybind11::module& m) {
  pybind11::class_<ParametrizedOptionListDescriptor, SettingDescriptor, std::shared_ptr<ParametrizedOptionListDescriptor>> parametrizedOptionListDescriptor(
      m, "ParametrizedOptionListDescriptor");

  parametrizedOptionListDescriptor.def(pybind11::init<std::string>(), pybind11::arg("description"));

  parametrizedOptionListDescriptor.def_property("default_value", &ParametrizedOptionListDescriptor::getDefaultOption,
                                                &ParametrizedOptionListDescriptor::setDefaultOption,
                                                "Default value of the settings");

  parametrizedOptionListDescriptor.def_property_readonly(
      "options", &ParametrizedOptionListDescriptor::getAllOptions,
      "Show all possible option strings and their matching settings");

  parametrizedOptionListDescriptor.def("add_option",
                                       pybind11::overload_cast<std::string>(&ParametrizedOptionListDescriptor::addOption),
                                       pybind11::arg("option"), "Add an option for which there are no specific settings");

  parametrizedOptionListDescriptor.def(
      "add_option", pybind11::overload_cast<std::string, DescriptorCollection>(&ParametrizedOptionListDescriptor::addOption),
      pybind11::arg("option"), pybind11::arg("settings"), "Add an option with specific attached settings");

  parametrizedOptionListDescriptor.def("option_settings", &ParametrizedOptionListDescriptor::getSettings,
                                       "Fetch settings corresponding to a particular option");

  parametrizedOptionListDescriptor.def_property_readonly("default_settings",
                                                         &ParametrizedOptionListDescriptor::getDefaultSettings,
                                                         "Fetch settings corresponding to the default option");
}

void init_int_list_descriptor(pybind11::module& m) {
  pybind11::class_<IntListDescriptor, SettingDescriptor, std::shared_ptr<IntListDescriptor>> intListDescriptor(
      m, "IntListDescriptor");

  intListDescriptor.def(pybind11::init<std::string>(), pybind11::arg("description"));

  intListDescriptor.def_property("default_value", &IntListDescriptor::getDefaultValue,
                                 &IntListDescriptor::setDefaultValue, "Default value of the settings");

  intListDescriptor.def_property("default_item_value", &IntListDescriptor::getDefaultItemValue,
                                 &IntListDescriptor::setDefaultItemValue, "Default item value for each item in the list");

  intListDescriptor.def_property("item_minimum", &IntListDescriptor::getItemMinimum, &IntListDescriptor::setItemMinimum,
                                 "Lower bound for items in the list");

  intListDescriptor.def_property("item_maximum", &IntListDescriptor::getItemMaximum, &IntListDescriptor::setItemMaximum,
                                 "Upper bound for items in the list");
}

void init_double_list_descriptor(pybind11::module& m) {
  pybind11::class_<DoubleListDescriptor, SettingDescriptor, std::shared_ptr<DoubleListDescriptor>> intListDescriptor(
      m, "DoubleListDescriptor");

  intListDescriptor.def(pybind11::init<std::string>(), pybind11::arg("description"));

  intListDescriptor.def_property("default_value", &DoubleListDescriptor::getDefaultValue,
                                 &DoubleListDescriptor::setDefaultValue, "Default value of the settings");

  intListDescriptor.def_property("default_item_value", &DoubleListDescriptor::getDefaultItemValue,
                                 &DoubleListDescriptor::setDefaultItemValue, "Default item value for each item in the list");

  intListDescriptor.def_property("item_minimum", &DoubleListDescriptor::getItemMinimum,
                                 &DoubleListDescriptor::setItemMinimum, "Lower bound for items in the list");

  intListDescriptor.def_property("item_maximum", &DoubleListDescriptor::getItemMaximum,
                                 &DoubleListDescriptor::setItemMaximum, "Upper bound for items in the list");
}

void init_string_list_descriptor(pybind11::module& m) {
  pybind11::class_<StringListDescriptor, SettingDescriptor, std::shared_ptr<StringListDescriptor>> stringListDescriptor(
      m, "StringListDescriptor");

  stringListDescriptor.def(pybind11::init<std::string>(), pybind11::arg("description"));

  stringListDescriptor.def_property("default_value", &StringListDescriptor::getDefaultValue,
                                    &StringListDescriptor::setDefaultValue, "Default value for the setting");

  stringListDescriptor.def_property("default_item_value", &StringListDescriptor::getDefaultItemValue,
                                    &StringListDescriptor::setDefaultItemValue,
                                    "Default value for each item in the string list");
}

void init_collection_list_descriptor(pybind11::module& m) {
  pybind11::class_<CollectionListDescriptor, SettingDescriptor, std::shared_ptr<CollectionListDescriptor>> collectionListDescriptor(
      m, "CollectionListDescriptor");

  collectionListDescriptor.def(pybind11::init<std::string, DescriptorCollection>(), pybind11::arg("description"),
                               pybind11::arg("base"));

  collectionListDescriptor.def_property_readonly("descriptor_collection", &CollectionListDescriptor::getDescriptorCollection);
}

using DescriptorValueVariant = RewriteAsVariant<GenericDescriptorMeta::Types>::Type;

GenericDescriptor convert(const DescriptorValueVariant& v) {
  return boost::apply_visitor([](const auto& t) { return GenericDescriptor(t); }, v);
}

void setItem(DescriptorCollection& coll, const std::string& key, const DescriptorValueVariant& v) {
  auto findIter =
      std::find_if(std::begin(coll), std::end(coll), [&](const auto& iterPair) -> bool { return iterPair.first == key; });

  if (findIter == std::end(coll)) {
    coll.push_back(key, convert(v));
  }
  else {
    findIter->second = convert(v);
  }
}

template<typename T>
void init_descriptor_collection(T& descriptor_collection) {
  descriptor_collection.def(pybind11::init<std::string>());

  descriptor_collection.def(
      "valid_value", pybind11::overload_cast<const ValueCollection&>(&DescriptorCollection::validValue, pybind11::const_),
      pybind11::arg("v"),
      R"delim(
      Checks whether a value collection matches the configuration specified by
      the members of this class.
    )delim");

  descriptor_collection.def(
      "__getitem__",
      [](DescriptorCollection& coll, const std::string& key) -> GenericDescriptorMeta::RefVariant {
        return GenericDescriptorMeta::convert(coll.get(key));
      },
      pybind11::return_value_policy::reference_internal);

  descriptor_collection.def("__setitem__", &setItem);

  descriptor_collection.def("__contains__", &DescriptorCollection::exists);

  descriptor_collection.def("__len__", &DescriptorCollection::size);
}

} // namespace

void init_setting_descriptors(pybind11::module& m) {
  init_setting_descriptor(m);

  pybind11::class_<DescriptorCollection, SettingDescriptor, std::shared_ptr<DescriptorCollection>> descriptor_collection(
      m, "DescriptorCollection",
      R"delim(
      Type-erased map-like container with string keys and Descriptor values.
    )delim");

  init_fundamental_descriptors(m);
  init_string_descriptor(m);
  init_file_descriptor(m);
  init_directory_descriptor(m);
  init_option_list_descriptor(m);
  init_parametrized_option_list_descriptor(m);
  init_int_list_descriptor(m);
  init_double_list_descriptor(m);
  init_string_list_descriptor(m);
  init_collection_list_descriptor(m);
  init_descriptor_collection(descriptor_collection);
}
