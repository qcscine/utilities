/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal Headers */
#include "Utils/UniversalSettings/ParametrizedOptionListDescriptor.h"
#include "Utils/UniversalSettings/Exceptions.h"
#include "Utils/UniversalSettings/ParametrizedOptionValue.h"

namespace Scine {
namespace Utils {
namespace UniversalSettings {

ParametrizedOptionListDescriptor::ParametrizedOptionListDescriptor(std::string propertyDescription)
  : SettingDescriptor(std::move(propertyDescription)) {
}

int ParametrizedOptionListDescriptor::optionCount() const {
  return static_cast<int>(options_.size());
}

bool ParametrizedOptionListDescriptor::optionExists(const std::string& name) const {
  return getIndex(name) != -1;
}

int ParametrizedOptionListDescriptor::getDefaultIndex() const {
  if (optionCount() == 0) {
    throw EmptyOptionListException(getPropertyDescription());
  }
  return defaultIndex_;
}

const std::string& ParametrizedOptionListDescriptor::getDefaultOption() const {
  return options_[getDefaultIndex()].first;
}

const std::vector<ParametrizedOptionListDescriptor::OptionAndSettings>& ParametrizedOptionListDescriptor::getAllOptions() const {
  return options_;
}

void ParametrizedOptionListDescriptor::addOption(std::string option) {
  DescriptorCollection emptyDescriptorCollection("(no settings required)");
  addOption(std::move(option), std::move(emptyDescriptorCollection));
}

void ParametrizedOptionListDescriptor::addOption(std::string option, DescriptorCollection optionSettings) {
  options_.emplace_back(std::move(option), std::move(optionSettings));
}

void ParametrizedOptionListDescriptor::setDefaultOption(const std::string& def) {
  int index = getIndex(def);
  if (index == -1) {
    throw OptionDoesNotExistException(def, getPropertyDescription());
  }
  defaultIndex_ = index;
}

const DescriptorCollection& ParametrizedOptionListDescriptor::getSettings(const std::string& option) const {
  int index = getIndex(option);
  if (index == -1) {
    throw OptionDoesNotExistException(option, getPropertyDescription());
  }
  return options_[index].second;
}

const DescriptorCollection& ParametrizedOptionListDescriptor::getDefaultSettings() const {
  return options_[getDefaultIndex()].second;
}

int ParametrizedOptionListDescriptor::getIndex(const std::string& option) const {
  int index = -1;
  for (int i = 0; i < optionCount(); ++i) {
    if (options_[i].first == option) {
      index = i;
    }
  }
  return index;
}

std::unique_ptr<SettingDescriptor> ParametrizedOptionListDescriptor::clone() const {
  return std::make_unique<ParametrizedOptionListDescriptor>(*this);
}

GenericValue ParametrizedOptionListDescriptor::getDefaultGenericValue() const {
  ParametrizedOptionValue v{getDefaultOption(), createDefaultValueCollection(getDefaultSettings())};
  return GenericValue::fromOptionWithSettings(std::move(v));
}

bool ParametrizedOptionListDescriptor::validValue(const GenericValue& v) const {
  if (!v.isOptionWithSettings()) {
    return false;
  }

  auto optionWithSettings = v.toOptionWithSettings();
  if (!optionExists(optionWithSettings.selectedOption)) {
    return false;
  }

  auto settings = getSettings(optionWithSettings.selectedOption);
  return settings.validValue(optionWithSettings.optionSettings);
}

inline std::string ParametrizedOptionListDescriptor::explainInvalidValue(const GenericValue& v) const {
  assert(!validValue(v));
  if (!v.isOptionWithSettings()) {
    return "Generic value for parametrized option list setting '" + getPropertyDescription() +
           "' is not a parametrized option list!";
  }
  auto optionWithSettings = v.toOptionWithSettings();
  if (!optionExists(optionWithSettings.selectedOption)) {
    return "Value " + optionWithSettings.selectedOption + " for parametrized option list setting '" +
           getPropertyDescription() + "' does not exist as an option!";
  }
  auto settings = getSettings(optionWithSettings.selectedOption);
  return settings.explainInvalidValue(optionWithSettings.optionSettings);
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
