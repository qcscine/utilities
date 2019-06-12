/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_YAML_H_
#define UTILS_YAML_H_

#include <Utils/Settings.h>
#include <yaml-cpp/yaml.h>

namespace Scine {
namespace Utils {

/**
 * @brief A custom exception for YAML parsing errors.
 */
class YAMLParsingException : public std::runtime_error {
 public:
  YAMLParsingException(std::string message) : std::runtime_error(message.c_str()) {
  }
  const char* what() const throw() {
    return std::runtime_error::what();
  }
};

/**
 * @brief Parses data from a YAML::Node (yaml-cpp) into the value collection of a settings object.
 *
 * @param settings The settings.
 * @param node     The read yaml node.
 * @param allowSuperfluous It true, the function will ignore kys in the node that are not present
 *                         in the settings object.
 *                         If false, the function will throw an error if such a case is encountered.
 */
static void nodeToSettings(Utils::Settings& settings, YAML::Node node, bool allowSuperfluous = false) {
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

#endif // UTILS_YAML_H_
