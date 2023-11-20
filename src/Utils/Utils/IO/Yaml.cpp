/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Yaml.h"
#include "Utils/MSVCCompatibility.h"
#include "Utils/UniversalSettings/GenericValueVariant.h"
#include <Utils/UniversalSettings/ParametrizedOptionValue.h>
#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/algorithm/iteration.hpp>
#include <boost/fusion/include/iteration.hpp>
#include <boost/fusion/include/std_tuple.hpp>
#include <boost/optional.hpp>
#include <cmath>

namespace Scine {
namespace Utils {
namespace {

boost::optional<int> safeInt(const std::string& source) {
  auto previousErrno = errno;
  errno = 0;
  char* interpretEnd = nullptr;
  auto ret = std::strtol(source.c_str(), &interpretEnd, 10);
  bool hasDigits = (source.find('.') != std::string::npos);
  bool failure = (errno == ERANGE || interpretEnd != source.c_str() + source.size() ||
                  ret > std::numeric_limits<int>::max() || hasDigits);
  errno = previousErrno;

  if (failure) {
    return boost::none;
  }

  return static_cast<int>(ret);
}

boost::optional<double> safeDouble(const std::string& source) {
  auto previousErrno = errno;
  errno = 0;
  char* interpretEnd = nullptr;
  double ret = std::strtod(source.c_str(), &interpretEnd);
  bool failure = (errno == ERANGE || interpretEnd != source.c_str() + source.size());
  errno = previousErrno;

  if (failure) {
    return boost::none;
  }

  return ret;
}

bool mapIsAParametrizedOption(const YAML::Node& node) {
  assert(node.IsMap());

  using namespace std::string_literals;
  static const std::vector<std::string> optionWithSettingsKeys = []() {
    std::vector<std::string> keys;
    keys.emplace_back("option_settings");
    keys.emplace_back("selected_option");
    return keys;
  }();

  std::vector<std::string> keys;
  for (auto subnode : node) {
    keys.push_back(subnode.first.as<std::string>());
  }
  std::sort(std::begin(keys), std::end(keys));

  assert(std::is_sorted(std::begin(optionWithSettingsKeys), std::end(optionWithSettingsKeys)));
  return keys == optionWithSettingsKeys;
}

void serialize(YAML::Emitter& emitter, const UniversalSettings::GenericValue& generic);

void serialize(YAML::Emitter& emitter, bool v) {
  emitter << v;
}

void serialize(YAML::Emitter& emitter, int v) {
  emitter << v;
}

void serialize(YAML::Emitter& emitter, const std::string& v) {
  emitter << v;
}

void serialize(YAML::Emitter& emitter, double v) {
  char str[80];
  double temp;
  if (modf(v, &temp) != 0) {
    sprintf(str, "%g", v);
    emitter << str;
  }
  else {
    sprintf(str, "%g.0", v);
    emitter << str;
  }
}

void serialize(YAML::Emitter& emitter, const UniversalSettings::ValueCollection& collection) {
  emitter << YAML::BeginMap;
  for (const auto& keyValuePair : collection) {
    emitter << YAML::Key << keyValuePair.first << YAML::Value;
    serialize(emitter, keyValuePair.second);
  }
  emitter << YAML::EndMap;
}

void serialize(YAML::Emitter& emitter, const UniversalSettings::ParametrizedOptionValue& option) {
  emitter << YAML::BeginMap;
  emitter << YAML::Key << "selected_option" << YAML::Value << option.selectedOption;
  emitter << YAML::Key << "option_settings" << YAML::Value;
  serialize(emitter, option.optionSettings);
  emitter << YAML::EndMap;
}

template<typename T>
void serialize(YAML::Emitter& emitter, const std::vector<T>& list) {
  emitter << YAML::BeginSeq;
  for (const auto& v : list) {
    serialize(emitter, v);
  }
  emitter << YAML::EndSeq;
}

void serialize(YAML::Emitter& emitter, const UniversalSettings::GenericValue& generic) {
  // clang-format off
  boost::fusion::for_each(
    Utils::UniversalSettings::GenericValueMeta::zip(
      Utils::UniversalSettings::GenericValueMeta::type_checkers(),
      Utils::UniversalSettings::GenericValueMeta::getters()
    ),
    [&](auto fns) {
      if (std::get<0>(fns)(generic)) {
        serialize(emitter, std::get<1>(fns)(generic));
      }
    }
  );
  // clang-format on
}

} // namespace

UniversalSettings::ValueCollection deserializeValueCollection(const YAML::Node& node) {
  UniversalSettings::ValueCollection collection;

  if (!node.IsMap()) {
    throw std::logic_error("Expected a YAML map node as root");
  }

  for (auto pair : node) {
    auto key = pair.first.as<std::string>();

    if (pair.second.IsScalar()) {
      // int / / double / bool / string
      const std::string& scalar = pair.second.Scalar();
      if (auto integer = safeInt(scalar)) {
        collection.addInt(key, integer.value());
      }
      else if (auto floating = safeDouble(scalar)) {
        collection.addDouble(key, floating.value());
      }
      else if (scalar == "true" || scalar == "True" || scalar == "TRUE") {
        collection.addBool(key, true);
      }
      else if (scalar == "false" || scalar == "False" || scalar == "FALSE") {
        collection.addBool(key, false);
      }
      else {
        collection.addString(key, scalar);
      }
    }
    else if (pair.second.IsSequence() && pair.second.size() > 0) {
      // List members
      YAML::Node firstElement = pair.second[0];
      if (firstElement.IsScalar()) {
        // Double-, IntList or StringList
        if (safeInt(firstElement.Scalar()) or safeDouble(firstElement.Scalar())) {
          bool isDoubleList = false;
          // Run one pass to check if Int- or DoubleList
          for (auto subnode : pair.second) {
            if (!(safeInt(subnode.Scalar()))) {
              isDoubleList = true;
              break;
            }
          }

          // Extract
          if (!isDoubleList) {
            std::vector<int> integers;
            for (auto subnode : pair.second) {
              integers.push_back(safeInt(subnode.Scalar()).value());
            }
            collection.addIntList(key, std::move(integers));
          }
          else {
            std::vector<double> doubles;
            for (auto subnode : pair.second) {
              doubles.push_back(safeDouble(subnode.Scalar()).value());
            }
            collection.addDoubleList(key, std::move(doubles));
          }
        }
        else {
          std::vector<std::string> strings;
          for (auto subnode : pair.second) {
            strings.push_back(subnode.Scalar());
          }
          collection.addStringList(key, std::move(strings));
        }
      }
      else {
        // Collection List
        std::vector<UniversalSettings::ValueCollection> collectionList;
        for (auto subnode : pair.second) {
          collectionList.push_back(deserializeValueCollection(subnode));
        }
        collection.addCollectionList(key, collectionList);
      }
    }
    else if (pair.second.IsMap()) {
      if (pair.second.size() == 2 && mapIsAParametrizedOption(pair.second)) {
        std::string selectedOption;
        Utils::UniversalSettings::ValueCollection optionSettings;

        for (auto subnodePair : pair.second) {
          const std::string key = subnodePair.first.as<std::string>();
          if (key == "selected_option") {
            if (!subnodePair.second.IsScalar()) {
              throw YAMLParsingException("Non-scalar type for expected string");
            }
            selectedOption = subnodePair.second.as<std::string>();
          }
          else {
            optionSettings = deserializeValueCollection(subnodePair.second);
          }
        }
        collection.addOptionWithSettings(
            key, Utils::UniversalSettings::ParametrizedOptionValue{std::move(selectedOption), std::move(optionSettings)});
      }
      else {
        collection.addCollection(key, deserializeValueCollection(pair.second));
      }
    }
  }

  return collection;
}

std::string yamlSerialize(const Utils::UniversalSettings::ValueCollection& collection) {
  YAML::Emitter emitter;
  emitter.SetMapFormat(YAML::Flow);
  emitter.SetSeqFormat(YAML::Flow);

  serialize(emitter, collection);
  return std::string{emitter.c_str()};
}

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
      else if (value.isDoubleList()) {
        auto val = it->second.as<std::vector<double>>();
        settings.modifyDoubleList(key, val);
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

void checkYamlKeyRecognition(const YAML::Node& node, std::vector<std::string> allowedKeywords) {
  for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
    auto key = it->first.as<std::string>();
    if (std::find(std::begin(allowedKeywords), std::end(allowedKeywords), key) == std::end(allowedKeywords)) {
      throw YAMLParsingException("Error: the key: '" + key + "' was not recognized.");
    }
  }
}

} // namespace Utils
} // namespace Scine
