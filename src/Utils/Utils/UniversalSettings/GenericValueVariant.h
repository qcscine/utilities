/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INCLUDE_UTILS_UNIVERSAL_SETTINGS_GENERIC_VALUE_VARIANT_H
#define INCLUDE_UTILS_UNIVERSAL_SETTINGS_GENERIC_VALUE_VARIANT_H

#include <Utils/UniversalSettings/GenericValue.h>
#include <Utils/UniversalSettings/ParametrizedOptionValue.h>
#include <boost/optional.hpp>
#include <boost/variant.hpp>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

template<typename Tuple>
struct RewriteAsVariant {
  template<std::size_t... Inds>
  static constexpr auto impl(std::index_sequence<Inds...> inds) -> boost::variant<std::tuple_element_t<Inds, Tuple>...>;

  using Type = decltype(impl(std::make_index_sequence<std::tuple_size<Tuple>::value>()));
};

/*! @brief Metaprogramming helper for GenericValue
 *
 * Groups types, factories, checkers and getter functions together for easier
 * reflection on all possible member types.
 *
 * Also includes several more general metaprogramming helper functions.
 *
 * Allows conversion of GenericValue to and from a variant of all its member
 * types. Very helpful for python bindings.
 *
 * Have a look at this file's impl for two examples of how to use this.
 */
struct GenericValueMeta {
  //! Member types of GenericValue
  using Types = std::tuple<bool, int, double, std::string, ValueCollection, ParametrizedOptionValue, GenericValue::IntList,
                           GenericValue::DoubleList, GenericValue::StringList, GenericValue::CollectionList>;

  //! GenericValue's equivalent as a variant (sum type)
  using Variant = RewriteAsVariant<Types>::Type;

  //! GenericValue's static factory functions from its member types
  static auto factories() {
    return std::make_tuple(&GenericValue::fromBool, &GenericValue::fromInt, &GenericValue::fromDouble,
                           &GenericValue::fromString, &GenericValue::fromCollection, &GenericValue::fromOptionWithSettings,
                           &GenericValue::fromIntList, &GenericValue::fromDoubleList, &GenericValue::fromStringList,
                           &GenericValue::fromCollectionList);
  }

  //! GenericValue's checker functions to see whether a wrapped type is a particular member
  static auto type_checkers() {
    return std::make_tuple(std::mem_fn(&GenericValue::isBool), std::mem_fn(&GenericValue::isInt),
                           std::mem_fn(&GenericValue::isDouble), std::mem_fn(&GenericValue::isString),
                           std::mem_fn(&GenericValue::isCollection), std::mem_fn(&GenericValue::isOptionWithSettings),
                           std::mem_fn(&GenericValue::isIntList), std::mem_fn(&GenericValue::isDoubleList),
                           std::mem_fn(&GenericValue::isStringList), std::mem_fn(&GenericValue::isCollectionList));
  }

  //! GenericValue's get functions retrieving a particular type
  static auto getters() {
    return std::make_tuple(std::mem_fn(&GenericValue::toBool), std::mem_fn(&GenericValue::toInt),
                           std::mem_fn(&GenericValue::toDouble), std::mem_fn(&GenericValue::toString),
                           std::mem_fn(&GenericValue::toCollection), std::mem_fn(&GenericValue::toOptionWithSettings),
                           std::mem_fn(&GenericValue::toIntList), std::mem_fn(&GenericValue::toDoubleList),
                           std::mem_fn(&GenericValue::toStringList), std::mem_fn(&GenericValue::toCollectionList));
  }

  template<std::size_t... Inds>
  static auto tuple_index_sequence_helper(std::index_sequence<Inds...> /*inds */) {
    return std::make_tuple(std::integral_constant<std::size_t, Inds>()...);
  }

  template<std::size_t N>
  static auto tuple_index_sequence() {
    return tuple_index_sequence_helper(std::make_index_sequence<N>());
  }

  template<class T, class U>
  static constexpr std::common_type_t<T, U> variadic_min(const T& t, const U& u) {
    return t < u ? t : u;
  }

  template<class First, class Second, class... Rest>
  static constexpr std::common_type_t<First, Second, Rest...> variadic_min(const First& f, const Second& s, const Rest&... t) {
    using ReturnType = std::common_type_t<First, Second, Rest...>;
    return variadic_min(variadic_min(static_cast<ReturnType>(f), static_cast<ReturnType>(s)), static_cast<ReturnType>(t)...);
  }

  template<std::size_t I, typename Tuple, std::size_t... Inds>
  static auto zipRow(const Tuple& tup, std::index_sequence<Inds...> /* inds */) {
    return std::make_tuple(std::get<I>(std::get<Inds>(tup))...);
  }

  template<typename Tuple, std::size_t... Inds>
  static auto zipHelper(const Tuple& tup, std::index_sequence<Inds...> /* inds */) {
    return std::make_tuple(zipRow<Inds>(tup, std::make_index_sequence<std::tuple_size<Tuple>::value>())...);
  }

  /*! @brief Variadic zipper for multiple tuples
   *
   * Rewrites the parameter pack of tuples into rows.
   *
   * Look at this file's impl for an example of how to use this.
   */
  template<typename... Tuples>
  static auto zip(Tuples&&... tuples) {
    constexpr std::size_t smallest = variadic_min(std::tuple_size<Tuples>::value...);
    return zipHelper(std::make_tuple(tuples...), std::make_index_sequence<smallest>());
  }

  //! Convert a generic value into a variant representation
  static Variant convert(const GenericValue& v);
  //! Convert a variant representation into a generic value
  static GenericValue convert(const Variant& v);

  //! Determine whether generic values are the same type or not
  static bool sameType(const GenericValue& a, const GenericValue& b);
};

} // namespace UniversalSettings
} // namespace Utils
} // namespace Scine

#endif
