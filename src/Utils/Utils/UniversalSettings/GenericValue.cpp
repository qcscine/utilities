/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal Headers */
#include "Utils/UniversalSettings/GenericValue.h"
#include "Utils/UniversalSettings/Exceptions.h"
#include "Utils/UniversalSettings/ParametrizedOptionValue.h"
#include "Utils/UniversalSettings/ValueCollection.h"
/* External Headers */
#include <boost/any.hpp>
#include <iostream>
namespace Scine {
namespace Utils {
namespace UniversalSettings {

struct GenericValue::Impl {
  boost::any value;
};

GenericValue::GenericValue() {
  pimpl_ = std::make_unique<Impl>();
}

GenericValue::~GenericValue() = default;

GenericValue::GenericValue(GenericValue&& /*rhs*/) noexcept = default;

GenericValue& GenericValue::operator=(GenericValue&& /*rhs*/) noexcept = default;

GenericValue::GenericValue(const GenericValue& rhs) : GenericValue() {
  pimpl_->value = rhs.pimpl_->value;
}

GenericValue& GenericValue::operator=(const GenericValue& rhs) {
  pimpl_->value = rhs.pimpl_->value;
  return *this;
}

GenericValue GenericValue::fromBool(bool v) {
  GenericValue gv;
  gv.pimpl_->value = v;
  return gv;
}

GenericValue GenericValue::fromInt(int v) {
  GenericValue gv;
  gv.pimpl_->value = v;
  return gv;
}

GenericValue GenericValue::fromIntList(IntList v) {
  GenericValue gv;
  gv.pimpl_->value = std::move(v);
  return gv;
}

GenericValue GenericValue::fromDouble(double v) {
  GenericValue gv;
  gv.pimpl_->value = v;
  return gv;
}

GenericValue GenericValue::fromString(std::string s) {
  GenericValue gv;
  gv.pimpl_->value = std::move(s);
  return gv;
}

GenericValue GenericValue::fromStringList(const GenericValue::StringList& v) {
  GenericValue gv;
  gv.pimpl_->value = v;
  return gv;
}

GenericValue GenericValue::fromCollection(ValueCollection s) {
  GenericValue gv;
  gv.pimpl_->value = std::move(s);
  return gv;
}

GenericValue GenericValue::fromCollectionList(GenericValue::CollectionList v) {
  GenericValue gv;
  gv.pimpl_->value = std::move(v);
  return gv;
}

GenericValue GenericValue::fromOptionWithSettings(ParametrizedOptionValue v) {
  GenericValue gv;
  gv.pimpl_->value = std::move(v);
  return gv;
}

bool GenericValue::isBool() const {
  return pimpl_->value.type() == typeid(bool);
}

bool GenericValue::isInt() const {
  return pimpl_->value.type() == typeid(int);
}

bool GenericValue::isIntList() const {
  return pimpl_->value.type() == typeid(IntList);
}

bool GenericValue::isDouble() const {
  return pimpl_->value.type() == typeid(double);
}

bool GenericValue::isString() const {
  // Clang 6.0.0 needed this ugly workaround.
  return (std::string(pimpl_->value.type().name()) == std::string(typeid(std::string).name()));
}

bool GenericValue::isStringList() const {
  return pimpl_->value.type() == typeid(StringList);
}

bool GenericValue::isCollection() const {
  return pimpl_->value.type() == typeid(ValueCollection);
}

bool GenericValue::isOptionWithSettings() const {
  return pimpl_->value.type() == typeid(ParametrizedOptionValue);
}

bool GenericValue::isCollectionList() const {
  return pimpl_->value.type() == typeid(CollectionList);
}

bool GenericValue::toBool() const {
  if (!isBool()) {
    throw InvalidValueConversionException{};
  }
  return boost::any_cast<bool>(pimpl_->value);
}

int GenericValue::toInt() const {
  if (!isInt()) {
    throw InvalidValueConversionException{};
  }
  return boost::any_cast<int>(pimpl_->value);
}

GenericValue::IntList GenericValue::toIntList() const {
  if (!isIntList()) {
    throw InvalidValueConversionException{};
  }
  return boost::any_cast<IntList>(pimpl_->value);
}

double GenericValue::toDouble() const {
  if (!isDouble()) {
    throw InvalidValueConversionException{};
  }
  return boost::any_cast<double>(pimpl_->value);
}

std::string GenericValue::toString() const {
  if (!isString()) {
    throw InvalidValueConversionException{};
  }
  return boost::any_cast<std::string>(pimpl_->value);
}

GenericValue::StringList GenericValue::toStringList() const {
  if (!isStringList()) {
    throw InvalidValueConversionException{};
  }
  return boost::any_cast<StringList>(pimpl_->value);
}

ValueCollection GenericValue::toCollection() const {
  if (!isCollection()) {
    throw InvalidValueConversionException{};
  }
  return boost::any_cast<ValueCollection>(pimpl_->value);
}

GenericValue::CollectionList GenericValue::toCollectionList() const {
  if (!isCollectionList()) {
    throw InvalidValueConversionException{};
  }
  return boost::any_cast<CollectionList>(pimpl_->value);
}

ParametrizedOptionValue GenericValue::toOptionWithSettings() const {
  if (!isOptionWithSettings()) {
    throw InvalidValueConversionException{};
  }
  return boost::any_cast<ParametrizedOptionValue>(pimpl_->value);
}
} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
