/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal Headers */
#include "Utils/UniversalSettings/InformationOutput.h"
#include "Utils/UniversalSettings/CollectionListDescriptor.h"
#include "Utils/UniversalSettings/DescriptorCollection.h"
#include "Utils/UniversalSettings/ParametrizedOptionListDescriptor.h"
/* External Headers */
#include <sstream>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

namespace {

/*! add quotes around text */
std::string quote(const std::string& text) {
  constexpr char quoteCharacter = '\"';
  return quoteCharacter + text + quoteCharacter;
}

} // namespace

void InformationOutput::print(const std::string& key, const DescriptorCollection& settings, std::ostream& out,
                              int indentation, bool outputCollectionTitle) {
  std::string indentationString;
  indentationString.append(indentation, ' ');

  if (outputCollectionTitle) {
    out << indentationString << quote(key) << " - " << quote(settings.getPropertyDescription())
        << ": settings collection" << std::endl;
  }

  for (const auto& d : settings) {
    const auto& settingKey = d.first;
    const auto& setting = d.second;
    auto type = setting.getType();
    out << indentationString << "- " << quote(settingKey) << " - " << quote(setting.getPropertyDescription()) << ": ";
    switch (type) {
      case GenericDescriptor::Type::Bool: {
        const auto& descriptor = setting.getBoolDescriptor();
        out << "boolean value. Default: " << std::boolalpha << descriptor.getDefaultValue() << std::endl;
        break;
      }
      case GenericDescriptor::Type::Int: {
        const auto& descriptor = setting.getIntDescriptor();
        out << "integer value. Bounds: [" << descriptor.getMinimum() << " - " << descriptor.getMaximum()
            << "]. Default: " << descriptor.getDefaultValue() << std::endl;
        break;
      }
      case GenericDescriptor::Type::Double: {
        const auto& descriptor = setting.getDoubleDescriptor();
        out << "floating-point value. Bounds: [" << descriptor.getMinimum() << " - " << descriptor.getMaximum()
            << "]. Default: " << descriptor.getDefaultValue() << std::endl;
        break;
      }
      case GenericDescriptor::Type::String: {
        const auto& descriptor = setting.getStringDescriptor();
        out << "string value. Default: " << quote(descriptor.getDefaultValue()) << std::endl;
        break;
      }
      case GenericDescriptor::Type::File: {
        const auto& descriptor = setting.getFileDescriptor();
        out << "file path. Default: " << quote(descriptor.getDefaultValue()) << std::endl;
        break;
      }
      case GenericDescriptor::Type::Directory: {
        const auto& descriptor = setting.getDirectoryDescriptor();
        out << "directory path. Default: " << quote(descriptor.getDefaultValue()) << std::endl;
        break;
      }
      case GenericDescriptor::Type::OptionList: {
        const auto& descriptor = setting.getOptionListDescriptor();
        out << "option list, with following possibilities: "
            << "(default: " << quote(descriptor.getDefaultValue()) << ")" << std::endl;
        for (const auto& option : descriptor.getAllOptions()) {
          out << indentationString << "    - " << quote(option) << std::endl;
        }
        break;
      }
      case GenericDescriptor::Type::SettingCollection: {
        const auto& descriptor = setting.getSettingCollectionDescriptor();
        out << "setting collection: " << std::endl;
        print(settingKey, descriptor, out, indentation + 4, false);
        break;
      }
      case GenericDescriptor::Type::ParametrizedOptionList: {
        const auto& descriptor = setting.getParametrizedOptionListDescriptor();
        out << "parametrized option list, with following possibilities: "
            << "(default: " << quote(descriptor.getDefaultOption()) << ")" << std::endl;
        for (const auto& option : descriptor.getAllOptions()) {
          out << indentationString << "    - " << quote(option.first) << ", with the following settings: " << std::endl;
          print("DUMMY", descriptor.getSettings(option.first), out, indentation + 8, false);
        }
        break;
      }
      case GenericDescriptor::Type::IntList: {
        const auto& descriptor = setting.getIntListDescriptor();
        std::stringstream items;
        for (const auto& i : descriptor.getDefaultValue()) {
          items << i << ", ";
        }
        out << "list of int values. "
            << "Default value for list: [" << items.str() << "], "
            << "bounds for element: [" << descriptor.getItemMinimum() << " - " << descriptor.getItemMaximum() << "], "
            << "default value for element: " << descriptor.getDefaultItemValue() << std::endl;
        break;
      }
      case GenericDescriptor::Type::DoubleList: {
        const auto& descriptor = setting.getDoubleListDescriptor();
        std::stringstream items;
        for (const auto& i : descriptor.getDefaultValue()) {
          items << i << ", ";
        }
        out << "list of double values. "
            << "Default value for list: [" << items.str() << "], "
            << "bounds for element: [" << descriptor.getItemMinimum() << " - " << descriptor.getItemMaximum() << "], "
            << "default value for element: " << descriptor.getDefaultItemValue() << std::endl;
        break;
      }
      case GenericDescriptor::Type::StringList: {
        const auto& descriptor = setting.getStringListDescriptor();
        std::stringstream items;
        for (const auto& i : descriptor.getDefaultValue()) {
          items << i << ", ";
        }
        out << "list of string values. "
            << "Default value for list: [" << items.str() << "], "
            << "default value for element: " << descriptor.getDefaultItemValue() << std::endl;
        break;
      }
      case GenericDescriptor::Type::CollectionList: {
        const auto& descriptor = setting.getCollectionListDescriptor();
        out << "list of collections: " << std::endl;
        print(settingKey, descriptor.getDescriptorCollection(), out, indentation + 4, false);
        break;
      }
      default: {
        throw Exception("GenericDescriptor not implemented.");
      }
    }
  }
}

void InformationOutput::printLong(const std::string& key, const DescriptorCollection& settings, std::ostream& out,
                                  int indentation) {
  std::string indentationString;
  indentationString.append(indentation, ' ');

  out << indentationString << "Setting collection, with key \"" << key << "\" and description \""
      << settings.getPropertyDescription() << "\"" << std::endl;

  for (const auto& d : settings) {
    const auto& settingKey = d.first;
    const auto& setting = d.second;
    auto type = setting.getType();
    out << indentationString << "- \"" << settingKey << "\": \"" << setting.getPropertyDescription() << "\"" << std::endl;
    out << indentationString << "  ";
    switch (type) {
      case GenericDescriptor::Type::Bool: {
        const auto& descriptor = setting.getBoolDescriptor();
        out << "Boolean value. Default: " << std::boolalpha << descriptor.getDefaultValue() << std::endl;
        break;
      }
      case GenericDescriptor::Type::Int: {
        const auto& descriptor = setting.getIntDescriptor();
        out << "Integer value. Bounds: between " << descriptor.getMinimum() << " and " << descriptor.getMaximum()
            << ". Default: " << descriptor.getDefaultValue() << std::endl;
        break;
      }
      case GenericDescriptor::Type::Double: {
        const auto& descriptor = setting.getDoubleDescriptor();
        out << "Floating-point value. Bounds: between " << descriptor.getMinimum() << " and " << descriptor.getMaximum()
            << ". Default: " << descriptor.getDefaultValue() << std::endl;
        break;
      }
      case GenericDescriptor::Type::String: {
        const auto& descriptor = setting.getStringDescriptor();
        out << "String value. Default: \"" << descriptor.getDefaultValue() << "\"" << std::endl;
        break;
      }
      case GenericDescriptor::Type::File: {
        const auto& descriptor = setting.getFileDescriptor();
        out << "File path. Default: \"" << descriptor.getDefaultValue() << "\"" << std::endl;
        break;
      }
      case GenericDescriptor::Type::Directory: {
        const auto& descriptor = setting.getDirectoryDescriptor();
        out << "Directory path. Default: \"" << descriptor.getDefaultValue() << "\"" << std::endl;
        break;
      }
      case GenericDescriptor::Type::OptionList: {
        const auto& descriptor = setting.getOptionListDescriptor();
        out << "Option list, with following possibilities: "
            << "(default: \"" << descriptor.getDefaultValue() << "\")" << std::endl;
        for (const auto& option : descriptor.getAllOptions()) {
          out << indentationString << "  - \"" << option << "\"" << std::endl;
        }
        break;
      }
      case GenericDescriptor::Type::SettingCollection: {
        const auto& descriptor = setting.getSettingCollectionDescriptor();
        out << "Setting collection with the following descriptors: " << std::endl;
        printLong(settingKey, descriptor, out, indentation + 2);
        break;
      }
      case GenericDescriptor::Type::ParametrizedOptionList: {
        const auto& descriptor = setting.getParametrizedOptionListDescriptor();
        out << "Parametrized option list, with following possibilities: "
            << "(default: \"" << descriptor.getDefaultOption() << "\")" << std::endl;
        for (const auto& option : descriptor.getAllOptions()) {
          out << indentationString << "  - \"" << option.first << "\", with the following settings: " << std::endl;
          printLong("DUMMY", descriptor.getSettings(option.first), out, indentation + 4);
        }
        break;
      }
      case GenericDescriptor::Type::IntList: {
        const auto& descriptor = setting.getIntListDescriptor();
        std::stringstream items;
        for (const auto& i : descriptor.getDefaultValue()) {
          items << i << ", ";
        }
        out << "list of int values. "
            << "Default value for list: [" << items.str() << "], "
            << "bounds for element: [" << descriptor.getItemMinimum() << " - " << descriptor.getItemMaximum() << "], "
            << "default value for element: " << descriptor.getDefaultItemValue() << std::endl;
        break;
      }
      case GenericDescriptor::Type::DoubleList: {
        const auto& descriptor = setting.getDoubleListDescriptor();
        std::stringstream items;
        for (const auto& i : descriptor.getDefaultValue()) {
          items << i << ", ";
        }
        out << "list of double values. "
            << "Default value for list: [" << items.str() << "], "
            << "bounds for element: [" << descriptor.getItemMinimum() << " - " << descriptor.getItemMaximum() << "], "
            << "default value for element: " << descriptor.getDefaultItemValue() << std::endl;
        break;
      }
      case GenericDescriptor::Type::StringList: {
        const auto& descriptor = setting.getStringListDescriptor();
        std::stringstream items;
        for (const auto& i : descriptor.getDefaultValue()) {
          items << i << ", ";
        }
        out << "list of string values. "
            << "Default value for list: [" << items.str() << "], "
            << "default value for element: " << descriptor.getDefaultItemValue() << std::endl;
        break;
      }
      case GenericDescriptor::Type::CollectionList: {
        const auto& descriptor = setting.getCollectionListDescriptor();
        out << "List of collections in which each item has the following descriptors: " << std::endl;
        print(settingKey, descriptor.getDescriptorCollection(), out, indentation + 4, false);
        break;
      }
      default: {
        throw Exception("GenericDescriptor not implemented.");
      }
    }
  }
}
} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
