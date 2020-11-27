/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal Headers */
#include "Utils/UniversalSettings/OptionListDescriptor.h"
#include "Utils/UniversalSettings/Exceptions.h"
#include "Utils/UniversalSettings/GenericValue.h"

namespace Scine {
namespace Utils {
namespace UniversalSettings {

OptionListDescriptor::OptionListDescriptor(std::string propertyDescription)
  : SettingDescriptor(std::move(propertyDescription)) {
}

int OptionListDescriptor::optionCount() const {
  return static_cast<int>(options_.size());
}

bool OptionListDescriptor::optionExists(const std::string& name) const {
  return getIndex(name) != -1;
}

int OptionListDescriptor::getDefaultIndex() const {
  if (optionCount() == 0) {
    throw EmptyOptionListException(getPropertyDescription());
  }
  return defaultIndex_;
}

const std::string& OptionListDescriptor::getDefaultOption() const {
  return options_[getDefaultIndex()];
}

void OptionListDescriptor::addOption(std::string option) {
  if (optionExists(option)) {
    throw OptionAlreadyExistsException(option, getPropertyDescription());
  }
  options_.push_back(std::move(option));
}

void OptionListDescriptor::setDefaultOption(const std::string& def) {
  int index = getIndex(def);
  if (index == -1) {
    throw OptionDoesNotExistException(def, getPropertyDescription());
  }
  defaultIndex_ = index;
}

int OptionListDescriptor::getIndex(const std::string& option) const {
  int index = -1;
  for (int i = 0; i < optionCount(); ++i) {
    if (options_[i] == option)
      index = i;
  }
  return index;
}

const std::vector<std::string>& OptionListDescriptor::getAllOptions() const {
  return options_;
}

std::unique_ptr<SettingDescriptor> OptionListDescriptor::clone() const {
  return std::make_unique<OptionListDescriptor>(*this);
}

GenericValue OptionListDescriptor::getDefaultGenericValue() const {
  return GenericValue::fromString(getDefaultOption());
}

bool OptionListDescriptor::validValue(const GenericValue& v) const {
  if (!v.isString()) {
    return false;
  }
  return optionExists(v.toString());
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
