/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal Headers */
#include "Utils/UniversalSettings/DescriptorCollection.h"
#include "Utils/UniversalSettings/Exceptions.h"
#include "Utils/UniversalSettings/ParametrizedOptionListDescriptor.h"
#include "Utils/UniversalSettings/ParametrizedOptionValue.h"
#include "Utils/UniversalSettings/ValueCollection.h"
/* External Headers */
#include <utility>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

DescriptorCollection::DescriptorCollection(std::string description) : SettingDescriptor(std::move(description)) {
}

DescriptorCollection::iterator DescriptorCollection::begin() {
  return descriptors_.begin();
}

DescriptorCollection::const_iterator DescriptorCollection::begin() const {
  return descriptors_.begin();
}

DescriptorCollection::iterator DescriptorCollection::end() {
  return descriptors_.end();
}

DescriptorCollection::const_iterator DescriptorCollection::end() const {
  return descriptors_.end();
}

void DescriptorCollection::push_back(std::string key, GenericDescriptor e) {
  for (const auto& d : descriptors_) {
    if (d.first == key) {
      throw AlreadyExistingDescriptorException(key);
    }
  }
  descriptors_.emplace_back(std::move(key), std::move(e));
}

GenericDescriptor& DescriptorCollection::operator[](unsigned i) {
  return descriptors_[i].second;
}

const GenericDescriptor& DescriptorCollection::operator[](unsigned i) const {
  return descriptors_[i].second;
}

GenericDescriptor& DescriptorCollection::at(unsigned i) {
  return descriptors_.at(i).second;
}

const GenericDescriptor& DescriptorCollection::at(unsigned i) const {
  return descriptors_.at(i).second;
}

bool DescriptorCollection::empty() const {
  return descriptors_.empty();
}

int DescriptorCollection::size() const {
  return static_cast<int>(descriptors_.size());
}

std::unique_ptr<SettingDescriptor> DescriptorCollection::clone() const {
  return std::make_unique<DescriptorCollection>(*this);
}

GenericValue DescriptorCollection::getDefaultGenericValue() const {
  return GenericValue::fromCollection(createDefaultValueCollection(*this));
}

GenericDescriptor& DescriptorCollection::get(const std::string& key) {
  return const_cast<GenericDescriptor&>(static_cast<const DescriptorCollection*>(this)->get(key));
}

const GenericDescriptor& DescriptorCollection::get(const std::string& key) const {
  for (const auto& d : descriptors_) {
    if (d.first == key) {
      return d.second;
    }
  }
  throw InexistingDescriptorException(key);
}

bool DescriptorCollection::exists(const std::string& key) const {
  for (const auto& descriptor : descriptors_) {
    if (descriptor.first == key) {
      return true;
    }
  }
  return false;
}

bool DescriptorCollection::validValue(const ValueCollection& v) const {
  // Check if every key in node is present in descriptors
  for (const auto& key : v.getKeys()) {
    if (!exists(key)) {
      return false;
    }
  }

  for (const auto& kvPair : descriptors_) {
    const auto& key = kvPair.first;

    if (!v.valueExists(key)) {
      return false;
    }

    const auto& genericDescriptor = kvPair.second;
    const auto& descriptor = genericDescriptor.getDescriptor();
    if (!descriptor.validValue(v.getValue(key))) {
      return false;
    }
  }

  return true;
}

bool DescriptorCollection::validValue(const GenericValue& v) const {
  if (!v.isCollection()) {
    return false;
  }

  return validValue(v.toCollection());
}

ValueCollection createDefaultValueCollection(const DescriptorCollection& descriptors) {
  ValueCollection defaultValues;
  for (const auto& d : descriptors) {
    const auto& key = d.first;
    const auto& setting = d.second;
    defaultValues.addGenericValue(key, setting.getDefaultValue());
  }
  return defaultValues;
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
