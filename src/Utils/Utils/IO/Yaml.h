/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_YAML_H_
#define UTILS_YAML_H_

#include <Utils/Settings.h>

namespace YAML {
class Node;
}
namespace Scine {
namespace Utils {

/**
 * @brief A custom exception for YAML parsing errors.
 */
class YAMLParsingException : public std::runtime_error {
 public:
  YAMLParsingException(const std::string& message) : std::runtime_error(message.c_str()) {
  }
  const char* what() const noexcept override {
    return std::runtime_error::what();
  }
};

/**
 * @brief Parses data from a YAML::Node (yaml-cpp) into a value collection
 */
Utils::UniversalSettings::ValueCollection deserializeValueCollection(const YAML::Node& node);

//! Serializes a value collection into YAML
std::string yamlSerialize(const Utils::UniversalSettings::ValueCollection& collection);

/**
 * @brief Parses data from a YAML::Node (yaml-cpp) into the value collection of a settings object.
 *
 * @param settings The settings.
 * @param node     The read yaml node.
 * @param allowSuperfluous It true, the function will ignore kys in the node that are not present
 *                         in the settings object.
 *                         If false, the function will throw an error if such a case is encountered.
 */
void nodeToSettings(Settings& settings, const YAML::Node& node, bool allowSuperfluous = false);

/**
 * @brief Checks whether all top level keys are recognized keywords.
 *
 * @param node The yaml node.
 * @param allowedKeywords The list of recognized keywords.
 * @throws YAMLParsingException if a key word is not recognized.
 */
void checkYamlKeyRecognition(const YAML::Node& node, std::vector<std::string> allowedKeywords);

} // namespace Utils
} // namespace Scine

#endif // UTILS_YAML_H_
