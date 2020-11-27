/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Yaml.h"
#include <yaml-cpp/yaml.h>

namespace Scine {
namespace Utils {

void nodeToSettings(Settings& settings, const YAML::Node& node, bool allowSuperfluous) {
  for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
    auto key = it->first.as<std::string>();
    if (settings.valueExists(key)) {
      auto value = settings.getValue(key);
      if (value.isInt()) {
        auto val = it->second.as<int>();
        settings.modifyInt(key, val);
      }
      else if (value.isBool()) {
        auto val = it->second.as<bool>();
        settings.modifyBool(key, val);
      }
      else if (value.isDouble()) {
        auto val = it->second.as<double>();
        settings.modifyDouble(key, val);
      }
      else if (value.isString()) {
        auto val = it->second.as<std::string>();
        settings.modifyString(key, val);
      }
      else if (value.isIntList()) {
        auto val = it->second.as<std::vector<int>>();
        settings.modifyIntList(key, val);
      }
      else if (value.isStringList()) {
        auto val = it->second.as<std::vector<std::string>>();
        settings.modifyStringList(key, val);
      }
      else if (value.isCollectionList()) {
        throw YAMLParsingException("Error: Cannot parse CollectionLists yet.");
      }
      else if (value.isCollection()) {
        throw YAMLParsingException("Error: Cannot parse Collections yet.");
      }
      else if (value.isOptionWithSettings()) {
        throw YAMLParsingException("Error: Cannot parse OptionWithSettings yet.");
      }
    }
    else if (!allowSuperfluous) {
      throw YAMLParsingException("Error: the key: '" + key + "' was not recognized.");
    }
  }
}

} // namespace Utils
} // namespace Scine
