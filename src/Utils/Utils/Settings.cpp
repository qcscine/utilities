/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Settings.h"
#include "Utils/UniversalSettings/CollectionListDescriptor.h"
#include "Utils/UniversalSettings/OptionListDescriptor.h"
#include "Utils/UniversalSettings/ParametrizedOptionListDescriptor.h"
#include "Utils/UniversalSettings/ParametrizedOptionValue.h"

namespace Scine {
namespace Utils {
namespace {

void transformLowercase(std::string& unknownCase) {
  std::transform(std::begin(unknownCase), std::end(unknownCase), std::begin(unknownCase),
                 [](unsigned char u) { return std::tolower(u); });
}

std::string makeLowercase(const std::string& unknownCase) {
  std::string lowercase;
  lowercase.resize(unknownCase.size());
  std::transform(std::begin(unknownCase), std::end(unknownCase), std::begin(lowercase),
                 [](unsigned char u) { return std::tolower(u); });
  return lowercase;
}

} // namespace

void Settings::normalizeStringCases(const UniversalSettings::DescriptorCollection& descriptors, ValueCollection& values) {
  using namespace UniversalSettings;

  for (const auto& keyDescriptorPair : descriptors) {
    const auto& key = keyDescriptorPair.first;
    if (!values.valueExists(key)) {
      // We're not validating, so this isn't an error!
      continue;
    }
    const auto& value = values.getValue(key);

    const auto& descriptor = keyDescriptorPair.second;
    const GenericDescriptor::Type type = descriptor.getType();
    switch (type) {
      case GenericDescriptor::Type::OptionList: {
        if (!value.isString()) {
          break;
        }

        std::string stringValue = value;
        const auto& optionList = descriptor.get<OptionListDescriptor>();
        if (optionList.optionExists(stringValue)) {
          break;
        }

        transformLowercase(stringValue);
        for (const auto& validOption : optionList.getAllOptions()) {
          const std::string lowercaseOption = makeLowercase(validOption);
          if (lowercaseOption == stringValue) {
            values.modifyValue(key, validOption);
            break;
          }
        }
        break;
      }
      case GenericDescriptor::Type::ParametrizedOptionList: {
        if (!value.isOptionWithSettings()) {
          break;
        }

        ParametrizedOptionValue optSettings = value.toOptionWithSettings();

        const auto& parametrizedOptionList = descriptor.get<ParametrizedOptionListDescriptor>();
        if (parametrizedOptionList.optionExists(optSettings.selectedOption)) {
          break;
        }

        transformLowercase(optSettings.selectedOption);
        for (const auto& validOptionDescriptorPair : parametrizedOptionList.getAllOptions()) {
          const std::string& validOptionIdentifier = validOptionDescriptorPair.first;
          const std::string lowercaseOption = makeLowercase(validOptionIdentifier);
          if (lowercaseOption == optSettings.selectedOption) {
            optSettings.selectedOption = validOptionIdentifier;
            // Recurse into DescriptorCollection + ValueCollection pair
            normalizeStringCases(parametrizedOptionList.getSettings(validOptionIdentifier), optSettings.optionSettings);
            values.modifyValue(key, optSettings);
            break;
          }
        }
        break;
      }
      case GenericDescriptor::Type::SettingCollection: {
        if (!value.isCollection()) {
          break;
        }

        ValueCollection subValues = value.toCollection();
        const auto& subDescriptors = descriptor.get<DescriptorCollection>();
        normalizeStringCases(subDescriptors, subValues);
        values.modifyValue(key, subValues);
        break;
      }
      case GenericDescriptor::Type::CollectionList: {
        if (!value.isCollectionList()) {
          break;
        }

        auto subValueCollections = value.toCollectionList();
        const auto& subDescriptors = descriptor.get<CollectionListDescriptor>();
        for (auto& subValueCollection : subValueCollections) {
          normalizeStringCases(subDescriptors.getDescriptorCollection(), subValueCollection);
        }
        values.modifyValue(key, subValueCollections);
        break;
      }
      default:
        break;
    }
  }
}

void Settings::normalizeStringCases() {
  normalizeStringCases(_fields, *this);
}

} // namespace Utils
} // namespace Scine
