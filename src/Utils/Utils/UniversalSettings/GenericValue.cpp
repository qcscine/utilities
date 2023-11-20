/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal Headers */
#include "Utils/UniversalSettings/Exceptions.h"
#include "Utils/UniversalSettings/GenericValueVariant.h"
#include "Utils/UniversalSettings/ParametrizedOptionValue.h"
#include "Utils/UniversalSettings/ValueCollection.h"
/* External Headers */
#include <boost/any.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/algorithm/iteration.hpp>
#include <boost/fusion/container/generation/make_map.hpp>
#include <boost/fusion/sequence/intrinsic/at_key.hpp>

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

GenericValue GenericValue::fromDoubleList(DoubleList v) {
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

bool GenericValue::isDoubleList() const {
  if (pimpl_->value.type() == typeid(DoubleList)) {
    return true;
  }
  // an empty list may still be a double list,
  // but boost assigns the type id of an int list
  return isEmptyIntList();
}

bool GenericValue::isEmptyIntList() const {
  if (pimpl_->value.type() == typeid(IntList)) {
    auto val = toIntList();
    return val.empty();
  }
  return false;
}

bool GenericValue::isDouble() const {
  return pimpl_->value.type() == typeid(double);
}

bool GenericValue::isString() const {
  // Clang 6.0.0 needed this ugly workaround.
  return (std::string(pimpl_->value.type().name()) == std::string(typeid(std::string).name()));
}

bool GenericValue::isStringList() const {
  if (pimpl_->value.type() == typeid(StringList)) {
    return true;
  }
  // an empty list may still be a string list,
  // but boost assigns the type id of an int list
  return isEmptyIntList();
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

GenericValue::DoubleList GenericValue::toDoubleList() const {
  if (!isDoubleList()) {
    throw InvalidValueConversionException{};
  }
  if (isEmptyIntList()) {
    return DoubleList{};
  }
  return boost::any_cast<DoubleList>(pimpl_->value);
}

double GenericValue::toDouble() const {
  if (!isDouble()) {
    throw InvalidValueConversionException{};
  }
  return boost::any_cast<double>(pimpl_->value);
}

std::string GenericValue::toString() const {
  if (isString()) {
    return boost::any_cast<std::string>(pimpl_->value);
  }
  if (isBool()) {
    auto v = boost::any_cast<bool>(pimpl_->value);
    return (v) ? "true" : "false";
  }
  if (isDouble()) {
    auto v = boost::any_cast<double>(pimpl_->value);
    return std::to_string(v);
  }
  if (isInt()) {
    auto v = boost::any_cast<int>(pimpl_->value);
    return std::to_string(v);
  }
  if (isIntList()) {
    auto v = boost::any_cast<IntList>(pimpl_->value);
    std::string ret = "[";
    for (const auto& i : v) {
      ret += std::to_string(i) + ", ";
    }
    ret = ret.substr(0, ret.size() - 2) + "]";
    return ret;
  }
  if (isDoubleList()) {
    auto v = boost::any_cast<DoubleList>(pimpl_->value);
    std::string ret = "[";
    for (const auto& i : v) {
      ret += std::to_string(i) + ", ";
    }
    ret = ret.substr(0, ret.size() - 2) + "]";
    return ret;
  }
  if (isStringList()) {
    auto v = boost::any_cast<StringList>(pimpl_->value);
    std::string ret = "[";
    for (const auto& i : v) {
      ret += i + ", ";
    }
    ret = ret.substr(0, ret.size() - 2) + "]";
    return ret;
  }
  if (isCollection()) {
    auto v = boost::any_cast<ValueCollection>(pimpl_->value);
    std::string ret = "{\n";
    for (const auto& [key, value] : v) {
      ret += "  " + key + ": " + value.toString() + ",\n";
    }
    ret = ret.substr(0, ret.size() - 2) + "\n}";
    return ret;
  }
  if (isCollectionList()) {
    std::string ret = "[";
    auto v = boost::any_cast<CollectionList>(pimpl_->value);
    for (const auto& coll : v) {
      std::string r = "{\n";
      for (const auto& [key, value] : coll) {
        r += "  " + key + ": " + value.toString() + ",\n";
      }
      r = r.substr(0, r.size() - 2) + "\n}";
      ret += r;
    }
    return ret;
  }
  throw InvalidValueConversionException{};
}

GenericValue::StringList GenericValue::toStringList() const {
  if (!isStringList()) {
    throw InvalidValueConversionException{};
  }
  if (isEmptyIntList()) {
    return StringList{};
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

bool operator==(const GenericValue& lhs, const GenericValue& rhs) {
  return boost::fusion::accumulate(GenericValueMeta::zip(GenericValueMeta::type_checkers(), GenericValueMeta::getters()),
                                   boost::optional<bool>(boost::none),
                                   [&](boost::optional<bool> maybeResult, auto t) -> boost::optional<bool> {
                                     if (!maybeResult) {
                                       auto&& check_fn = std::get<0>(t);
                                       auto&& get_fn = std::get<1>(t);
                                       if (check_fn(lhs)) {
                                         return check_fn(rhs) && get_fn(lhs) == get_fn(rhs);
                                       }
                                     }

                                     return maybeResult;
                                   })
      .value();
}

bool operator!=(const GenericValue& lhs, const GenericValue& rhs) {
  return !(lhs == rhs);
}

/* Implementation of GenericValue's templated member functions over its member
 * types
 */

namespace {

template<typename T, std::size_t... Inds>
auto makeMapImpl(const T& matchingTuple, std::index_sequence<Inds...> /* inds */) {
  return boost::fusion::make_map<std::tuple_element_t<Inds, GenericValueMeta::Types>...>(std::get<Inds>(matchingTuple)...);
}

template<typename T>
auto makeMap(const T& matchingTuple) {
  constexpr unsigned N = std::tuple_size<GenericValueMeta::Types>::value;
  return makeMapImpl(matchingTuple, std::make_index_sequence<N>());
}

} // namespace

static_assert(std::tuple_size<GenericValue::Types>::value == 10,
              "Have you changed the possible member types of GenericValue? You need to fix the template instantiations "
              "below, too");

// Constructor from members
template<typename T, std::enable_if_t<Detail::TypeInTuple<std::decay_t<T>, GenericValue::Types>::value, bool*>>
GenericValue::GenericValue(T t) {
  auto factory = boost::fusion::at_key<std::decay_t<T>>(makeMap(GenericValueMeta::factories()));
  *this = factory(std::move(t));
}

template GenericValue::GenericValue(bool);
template GenericValue::GenericValue(int);
template GenericValue::GenericValue(double);
template GenericValue::GenericValue(std::string);
template GenericValue::GenericValue(ValueCollection);
template GenericValue::GenericValue(ParametrizedOptionValue);
template GenericValue::GenericValue(IntList);
template GenericValue::GenericValue(DoubleList);
template GenericValue::GenericValue(StringList);
template GenericValue::GenericValue(CollectionList);

GenericValue::GenericValue(const char* str) : GenericValue(std::string{str}) {
}

// GenericValue assignment operator from members
template<typename T>
std::enable_if_t<Detail::TypeInTuple<std::decay_t<T>, GenericValue::Types>::value, GenericValue&>
GenericValue::operator=(T&& t) {
  auto factory = boost::fusion::at_key<std::decay_t<T>>(makeMap(GenericValueMeta::factories()));
  *this = factory(std::forward<T>(t));
  return *this;
}

template GenericValue& GenericValue::operator=(bool&&);
template GenericValue& GenericValue::operator=(int&&);
template GenericValue& GenericValue::operator=(double&&);
template GenericValue& GenericValue::operator=(std::string&&);
template GenericValue& GenericValue::operator=(ValueCollection&&);
template GenericValue& GenericValue::operator=(ParametrizedOptionValue&&);
template GenericValue& GenericValue::operator=(IntList&&);
template GenericValue& GenericValue::operator=(DoubleList&&);
template GenericValue& GenericValue::operator=(StringList&&);
template GenericValue& GenericValue::operator=(CollectionList&&);

GenericValue& GenericValue::operator=(const char* str) {
  *this = fromString(std::string{str});
  return *this;
}

// Equality checker with member types
template<typename T>
std::enable_if_t<Detail::TypeInTuple<std::decay_t<T>, GenericValue::Types>::value, bool>
GenericValue::operator==(const T& other) const {
  auto checker = boost::fusion::at_key<std::decay_t<T>>(makeMap(GenericValueMeta::type_checkers()));
  auto converter = boost::fusion::at_key<std::decay_t<T>>(makeMap(GenericValueMeta::getters()));
  return checker(*this) && converter(*this) == other;
}

template bool GenericValue::operator==(const bool&) const;
template bool GenericValue::operator==(const int&) const;
template bool GenericValue::operator==(const double&) const;
template bool GenericValue::operator==(const std::string&) const;
template bool GenericValue::operator==(const ValueCollection&) const;
template bool GenericValue::operator==(const ParametrizedOptionValue&) const;
template bool GenericValue::operator==(const IntList&) const;
template bool GenericValue::operator==(const DoubleList&) const;
template bool GenericValue::operator==(const StringList&) const;
template bool GenericValue::operator==(const CollectionList&) const;

bool GenericValue::operator==(const char* str) const {
  return isString() && toString() == str;
}

bool GenericValue::operator!=(const char* str) const {
  return !(*this == str);
}

// Implicit conversion to member type
template<typename T, std::enable_if_t<Detail::TypeInTuple<std::decay_t<T>, GenericValue::Types>::value, bool*>>
GenericValue::operator T() const {
  auto checker = boost::fusion::at_key<std::decay_t<T>>(makeMap(GenericValueMeta::type_checkers()));
  auto converter = boost::fusion::at_key<std::decay_t<T>>(makeMap(GenericValueMeta::getters()));
  if (!checker(*this)) {
    throw std::runtime_error("GenericValue is not the type being implicitly casted to!");
  }

  return converter(*this);
}

template GenericValue::operator bool() const;
template GenericValue::operator int() const;
template GenericValue::operator double() const;
template GenericValue::operator std::string() const;
template GenericValue::operator ValueCollection() const;
template GenericValue::operator ParametrizedOptionValue() const;
template GenericValue::operator IntList() const;
template GenericValue::operator DoubleList() const;
template GenericValue::operator StringList() const;
template GenericValue::operator CollectionList() const;

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
