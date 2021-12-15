/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal Headers */
#include "Utils/UniversalSettings/ValueCollection.h"
#include "Utils/UniversalSettings/Exceptions.h"
#include "Utils/UniversalSettings/ParametrizedOptionValue.h"
/* External Headers */
#include <algorithm>
#include <cassert>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

void ValueCollection::addGenericValue(std::string name, GenericValue value) {
  assert(!name.empty());
  if (valueExists(name)) {
    throw AlreadyExistingValueException(name);
  }
  values_.emplace_back(std::move(name), std::move(value));
}

GenericValue ValueCollection::getValue(const std::string& name) const {
  return getGenericValue(name);
}

void ValueCollection::modifyValue(const std::string& name, GenericValue value) {
  auto it = findName(name);
  if (it == values_.end()) {
    throw InexistingValueException(name);
  }
  it->second = std::move(value);
}

void ValueCollection::dropValue(const std::string& name) {
  auto it = findName(name);
  if (it == values_.end()) {
    throw InexistingValueException(name);
  }
  values_.erase(
      std::remove_if(std::begin(values_), std::end(values_), [&](const auto& pair) { return pair.first == name; }),
      std::end(values_));
}

void ValueCollection::addBool(std::string name, bool value) {
  addGenericValue(std::move(name), GenericValue::fromBool(value));
}

void ValueCollection::addInt(std::string name, int value) {
  addGenericValue(std::move(name), GenericValue::fromInt(value));
}

void ValueCollection::addDouble(std::string name, double value) {
  addGenericValue(std::move(name), GenericValue::fromDouble(value));
}

void ValueCollection::addString(std::string name, std::string value) {
  addGenericValue(std::move(name), GenericValue::fromString(std::move(value)));
}

void ValueCollection::addCollection(std::string name, ValueCollection value) {
  addGenericValue(std::move(name), GenericValue::fromCollection(std::move(value)));
}

void ValueCollection::addOptionWithSettings(std::string name, ParametrizedOptionValue value) {
  addGenericValue(std::move(name), GenericValue::fromOptionWithSettings(std::move(value)));
}

void ValueCollection::modifyBool(const std::string& name, bool value) {
  const auto& genericValue = getGenericValue(name);
  if (!genericValue.isBool()) {
    throw InvalidValueConversionException();
  }
  modifyValue(name, GenericValue::fromBool(value));
}

void ValueCollection::modifyInt(const std::string& name, int value) {
  const auto& genericValue = getGenericValue(name);
  if (!genericValue.isInt()) {
    throw InvalidValueConversionException();
  }
  modifyValue(name, GenericValue::fromInt(value));
}

void ValueCollection::modifyDouble(const std::string& name, double value) {
  const auto& genericValue = getGenericValue(name);
  if (!genericValue.isDouble()) {
    throw InvalidValueConversionException();
  }
  modifyValue(name, GenericValue::fromDouble(value));
}

void ValueCollection::modifyString(const std::string& name, std::string value) {
  const auto& genericValue = getGenericValue(name);
  if (!genericValue.isString()) {
    throw InvalidValueConversionException();
  }
  modifyValue(name, GenericValue::fromString(std::move(value)));
}

void ValueCollection::modifyCollection(const std::string& name, ValueCollection value) {
  const auto& genericValue = getGenericValue(name);
  if (!genericValue.isCollection()) {
    throw InvalidValueConversionException();
  }
  modifyValue(name, GenericValue::fromCollection(std::move(value)));
}

void ValueCollection::modifyOptionsWithSettings(const std::string& name, ParametrizedOptionValue value) {
  const auto& genericValue = getGenericValue(name);
  if (!genericValue.isOptionWithSettings()) {
    throw InvalidValueConversionException();
  }
  modifyValue(name, GenericValue::fromOptionWithSettings(std::move(value)));
}

bool ValueCollection::getBool(const std::string& name) const {
  const auto& genericValue = getGenericValue(name);
  try {
    return genericValue.toBool();
  }
  catch (InvalidValueConversionException&) {
    throw ValueHasDifferentTypeException(name);
  }
}

int ValueCollection::getInt(const std::string& name) const {
  const auto& genericValue = getGenericValue(name);
  try {
    return genericValue.toInt();
  }
  catch (InvalidValueConversionException&) {
    throw ValueHasDifferentTypeException(name);
  }
}

double ValueCollection::getDouble(const std::string& name) const {
  const auto& genericValue = getGenericValue(name);
  try {
    return genericValue.toDouble();
  }
  catch (InvalidValueConversionException&) {
    throw ValueHasDifferentTypeException(name);
  }
}

std::string ValueCollection::getString(const std::string& name) const {
  const auto& genericValue = getGenericValue(name);
  try {
    return genericValue.toString();
  }
  catch (InvalidValueConversionException&) {
    throw ValueHasDifferentTypeException(name);
  }
}

ValueCollection ValueCollection::getCollection(const std::string& name) const {
  const auto& genericValue = getGenericValue(name);
  try {
    return genericValue.toCollection();
  }
  catch (InvalidValueConversionException&) {
    throw ValueHasDifferentTypeException(name);
  }
}

ParametrizedOptionValue ValueCollection::getOptionWithSettings(const std::string& name) const {
  const auto& genericValue = getGenericValue(name);
  try {
    return genericValue.toOptionWithSettings();
  }
  catch (InvalidValueConversionException&) {
    throw ValueHasDifferentTypeException(name);
  }
}

void ValueCollection::addIntList(std::string name, GenericValue::IntList value) {
  addGenericValue(std::move(name), GenericValue::fromIntList(std::move(value)));
}

GenericValue::IntList ValueCollection::getIntList(const std::string& name) const {
  const auto& genericValue = getGenericValue(name);
  try {
    return genericValue.toIntList();
  }
  catch (InvalidValueConversionException&) {
    throw ValueHasDifferentTypeException(name);
  }
}

void ValueCollection::modifyIntList(const std::string& name, GenericValue::IntList value) {
  const auto& genericValue = getGenericValue(name);
  if (!genericValue.isIntList()) {
    throw InvalidValueConversionException();
  }
  modifyValue(name, GenericValue::fromIntList(std::move(value)));
}

void ValueCollection::addDoubleList(std::string name, GenericValue::DoubleList value) {
  addGenericValue(std::move(name), GenericValue::fromDoubleList(std::move(value)));
}

GenericValue::DoubleList ValueCollection::getDoubleList(const std::string& name) const {
  const auto& genericValue = getGenericValue(name);
  try {
    return genericValue.toDoubleList();
  }
  catch (InvalidValueConversionException&) {
    throw ValueHasDifferentTypeException(name);
  }
}

void ValueCollection::modifyDoubleList(const std::string& name, GenericValue::DoubleList value) {
  const auto& genericValue = getGenericValue(name);
  if (!genericValue.isDoubleList()) {
    throw InvalidValueConversionException();
  }
  modifyValue(name, GenericValue::fromDoubleList(std::move(value)));
}

void ValueCollection::addStringList(std::string name, GenericValue::StringList value) {
  addGenericValue(std::move(name), GenericValue::fromStringList(std::move(value)));
}

GenericValue::StringList ValueCollection::getStringList(const std::string& name) const {
  const auto& genericValue = getGenericValue(name);
  try {
    return genericValue.toStringList();
  }
  catch (InvalidValueConversionException&) {
    throw ValueHasDifferentTypeException(name);
  }
}

void ValueCollection::modifyStringList(const std::string& name, GenericValue::StringList value) {
  const auto& genericValue = getGenericValue(name);
  if (!genericValue.isStringList()) {
    throw InvalidValueConversionException();
  }
  modifyValue(name, GenericValue::fromStringList(std::move(value)));
}

void ValueCollection::addCollectionList(std::string name, GenericValue::CollectionList value) {
  addGenericValue(std::move(name), GenericValue::fromCollectionList(std::move(value)));
}

GenericValue::CollectionList ValueCollection::getCollectionList(const std::string& name) const {
  const auto& genericValue = getGenericValue(name);
  try {
    return genericValue.toCollectionList();
  }
  catch (InvalidValueConversionException&) {
    throw ValueHasDifferentTypeException(name);
  }
}

void ValueCollection::modifyCollectionList(const std::string& name, GenericValue::CollectionList value) {
  const auto& genericValue = getGenericValue(name);
  if (!genericValue.isCollectionList()) {
    throw InvalidValueConversionException();
  }
  modifyValue(name, GenericValue::fromCollectionList(std::move(value)));
}

ValueCollection::Container::const_iterator ValueCollection::findName(const std::string& name) const {
  auto it = std::find_if(values_.begin(), values_.end(), [&name](const auto& p) { return p.first == name; });
  return it;
}

ValueCollection::Container::iterator ValueCollection::findName(const std::string& name) {
  auto it = std::find_if(values_.begin(), values_.end(), [&name](const auto& p) { return p.first == name; });
  return it;
}

const GenericValue& ValueCollection::getGenericValue(const std::string& name) const {
  auto it = findName(name);
  if (it == values_.end())
    throw InexistingValueException(name);

  const auto& genericValue = it->second;
  return genericValue;
}

int ValueCollection::count() const {
  return static_cast<int>(values_.size());
}

int ValueCollection::size() const {
  return static_cast<int>(values_.size());
}

bool ValueCollection::empty() const {
  return values_.empty();
}

std::vector<std::string> ValueCollection::getKeys() const {
  std::vector<std::string> keys;
  for (const auto& v : values_) {
    keys.push_back(v.first);
  }
  return keys;
}

bool ValueCollection::valueExists(const std::string& name) const {
  return findName(name) != values_.end();
}

inline bool equalToDefault(const std::string& key, const GenericValue& value) {
  auto descriptor = OptionListDescriptor(key);
  return value == descriptor.getDefaultGenericValue();
}

bool operator==(const ValueCollection& lhs, const ValueCollection& rhs) {
  if (lhs.empty() && rhs.empty()) {
    return true;
  }
  auto lhsKeys = lhs.getKeys();
  auto rhsKeys = rhs.getKeys();
  for (const auto& key : lhsKeys) {
    // key not present or key present in both, but value is different
    if (!rhs.valueExists(key) || lhs.getValue(key) != rhs.getValue(key)) {
      return false;
    }
  }

  // identical the other way around
  for (const auto& key : rhsKeys) {
    if (!lhs.valueExists(key) || lhs.getValue(key) != rhs.getValue(key)) {
      return false;
    }
  }
  // all checks survived -> settings are equal
  return true;
}

bool operator!=(const ValueCollection& lhs, const ValueCollection& rhs) {
  return !(lhs == rhs);
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
