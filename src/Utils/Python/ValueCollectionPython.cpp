/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "boost/bind.hpp"
#include "boost/variant.hpp"
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/ParametrizedOptionValue.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/algorithm/iteration.hpp>
#include <boost/fusion/include/iteration.hpp>
#include <boost/fusion/include/std_tuple.hpp>
#include <boost/optional.hpp>

using namespace Scine::Utils;

namespace pybind11 {
namespace detail {

template<typename... Ts>
struct type_caster<boost::variant<Ts...>> : variant_caster<boost::variant<Ts...>> {};

template<>
struct visit_helper<boost::variant> {
  template<typename... Args>
  static auto call(Args&&... args) -> decltype(boost::apply_visitor(args...)) {
    return boost::apply_visitor(args...);
  }
};

} // namespace detail
} // namespace pybind11

namespace {

template<typename Tuple>
struct VariantFromTuple {
  template<std::size_t... Inds>
  static constexpr auto impl(std::index_sequence<Inds...> inds) -> boost::variant<std::tuple_element_t<Inds, Tuple>...>;

  using Type = decltype(impl(std::make_index_sequence<std::tuple_size<Tuple>::value>()));
};

using namespace UniversalSettings;
struct GenericValueMeta {
  using Types = std::tuple<bool, int, double, std::string, ValueCollection, ParametrizedOptionValue,
                           GenericValue::IntList, GenericValue::StringList, GenericValue::CollectionList>;

  using Variant = VariantFromTuple<Types>::Type;

  static auto factories() {
    return std::make_tuple(&GenericValue::fromBool, &GenericValue::fromInt, &GenericValue::fromDouble,
                           &GenericValue::fromString, &GenericValue::fromCollection, &GenericValue::fromOptionWithSettings,
                           &GenericValue::fromIntList, &GenericValue::fromStringList, &GenericValue::fromCollectionList);
  }

  static auto type_checkers() {
    return std::make_tuple(std::mem_fn(&GenericValue::isBool), std::mem_fn(&GenericValue::isInt),
                           std::mem_fn(&GenericValue::isDouble), std::mem_fn(&GenericValue::isString),
                           std::mem_fn(&GenericValue::isCollection), std::mem_fn(&GenericValue::isOptionWithSettings),
                           std::mem_fn(&GenericValue::isIntList), std::mem_fn(&GenericValue::isStringList),
                           std::mem_fn(&GenericValue::isCollectionList));
  }

  static auto getters() {
    return std::make_tuple(std::mem_fn(&GenericValue::toBool), std::mem_fn(&GenericValue::toInt),
                           std::mem_fn(&GenericValue::toDouble), std::mem_fn(&GenericValue::toString),
                           std::mem_fn(&GenericValue::toCollection), std::mem_fn(&GenericValue::toOptionWithSettings),
                           std::mem_fn(&GenericValue::toIntList), std::mem_fn(&GenericValue::toStringList),
                           std::mem_fn(&GenericValue::toCollectionList));
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

  template<typename... Tuples>
  static auto zip(Tuples&&... tuples) {
    constexpr std::size_t smallest = variadic_min(std::tuple_size<Tuples>::value...);
    return zipHelper(std::make_tuple(tuples...), std::make_index_sequence<smallest>());
  }

  static auto check_convert_pairs() {
    return zip(type_checkers(), getters());
  }

  static auto index_factory_pairs() {
    return zip(tuple_index_sequence<std::tuple_size<Types>::value>(), factories());
  }

  struct ToVariantConverter {
    using result_type = boost::optional<Variant>;

    template<typename CheckGetPair>
    result_type operator()(result_type maybeVariant, const CheckGetPair& t) const {
      if (maybeVariant) {
        return maybeVariant;
      }

      auto&& check_fn = std::get<0>(t);
      auto&& getter_fn = std::get<1>(t);
      if (check_fn(ref)) {
        return Variant(getter_fn(ref));
      }

      return maybeVariant;
    }

    ToVariantConverter(const GenericValue& f) : ref(f) {
    }
    const GenericValue& ref;
  };

  static Variant convert(const GenericValue& v) {
    auto maybeVariant = boost::fusion::accumulate(check_convert_pairs(), ToVariantConverter::result_type(boost::none),
                                                  ToVariantConverter(v));
    return maybeVariant.value();
  }

  struct ToGenericValueConverter {
    using result_type = boost::optional<GenericValue>;

    template<typename IndexFactoryPair>
    result_type operator()(result_type maybeGenericValue, const IndexFactoryPair& t) const {
      if (maybeGenericValue) {
        return maybeGenericValue;
      }

      const auto type_index_constant = std::get<0>(t);
      auto&& factory_fn = std::get<1>(t);
      constexpr std::size_t variantTypeIndex = decltype(type_index_constant)::value;
      if (ref.which() == variantTypeIndex) {
        using Subtype = std::tuple_element_t<variantTypeIndex, Types>;
        return factory_fn(boost::get<Subtype>(ref));
      }

      return maybeGenericValue;
    }

    ToGenericValueConverter(const Variant& v) : ref(v) {
    }
    const Variant& ref;
  };

  static GenericValue convert(const Variant& v) {
    auto maybeGenericValue = boost::fusion::accumulate(
        index_factory_pairs(), ToGenericValueConverter::result_type(boost::none), ToGenericValueConverter(v));
    return maybeGenericValue.value();
  }
};

struct Updater {
  using PybindTypes =
      std::tuple<pybind11::bool_, pybind11::int_, pybind11::float_, pybind11::str, ValueCollection,
                 ParametrizedOptionValue, GenericValue::IntList, GenericValue::StringList, GenericValue::CollectionList>;
  static_assert(std::tuple_size<PybindTypes>::value == std::tuple_size<GenericValueMeta::Types>::value,
                "Modified list of types in GenericValueMeta without adapting PybindTypes in update fn");

  template<typename T>
  void operator()(T&& t) const {
    auto type_index_constant = std::get<0>(t);
    auto&& factory_fn = std::get<1>(t);
    auto&& typecheck_fn = std::get<2>(t);

    constexpr unsigned typeIndex = decltype(type_index_constant)::value;
    using PyType = std::tuple_element_t<typeIndex, PybindTypes>;
    using CppType = std::tuple_element_t<typeIndex, GenericValueMeta::Types>;
    if (pybind11::isinstance<PyType>(handle)) {
      const auto dictValue = handle.cast<CppType>();
      GenericValue wrapped = factory_fn(dictValue);
      if (collection.valueExists(key)) {
        if (preserveTypes && !typecheck_fn(collection.getValue(key))) {
          throw std::runtime_error("GenericValue already exists as a different type!");
        }
        collection.modifyValue(key, wrapped);
      }
      else {
        collection.addGenericValue(key, wrapped);
      }
    }
  }

  Updater(const std::string& k, pybind11::handle item, ValueCollection& v, bool preserve)
    : key(k), handle(item), collection(v), preserveTypes(preserve) {
  }

  const std::string& key;
  pybind11::handle handle;
  ValueCollection& collection;
  bool preserveTypes;
};

pybind11::dict toDict(const ValueCollection& coll) {
  pybind11::dict dictionary;
  for (const auto& key : coll.getKeys()) {
    dictionary[key.c_str()] = GenericValueMeta::convert(coll.getValue(key));
  }
  return dictionary;
}

} // namespace

/* This function is outside the anonymous namespace so that we can reuse it for
 * a Settings constructor with a python dictionary
 */
void update(ValueCollection& coll, const pybind11::dict& dictionary, bool preserveTypes) {
  for (auto item : dictionary) {
    if (!pybind11::isinstance<pybind11::str>(item.first)) {
      throw std::runtime_error("Keys in the update dictionary must be strings!");
    }
    const auto key = std::string(pybind11::str(item.first));

    boost::fusion::for_each(
        GenericValueMeta::zip(GenericValueMeta::tuple_index_sequence<std::tuple_size<GenericValueMeta::Types>::value>(),
                              GenericValueMeta::factories(), GenericValueMeta::type_checkers()),
        Updater(key, item.second, coll, preserveTypes));
  }
}

void init_value_collection(pybind11::module& m) {
  pybind11::class_<ValueCollection, std::shared_ptr<ValueCollection>> value_collection(m, "ValueCollection",
                                                                                       R"delim(
      Type-erased C++ map-like container with string keys that can store the
      following types of values: ``bool``, ``int``, ``float``, ``str``,
      ``ValueCollection`` (enables nesting!), ``List[int]``, ``List[str]``,
      and ``List[ValueCollection]``.

      Has members to imitate behavior of a python dictionary with string keys.

      >>> coll = ValueCollection()
      >>> coll
      {}
      >>> coll["a"] = [1, 2, 3, 4]  # Add or lookup elements by key
      >>> coll
      {'a': [1, 2, 3, 4]}
      >>> "b" in coll
      False
      >>> coll.get("b", 4.0)  # Default values for missing keys
      4.0
      >>> list(coll)  # Iteration is supported, yielding key-value tuples
      [('a', [1, 2, 3, 4])]
      >>> coll.update({"b": 4.0})
      >>> len(coll)
      2
    )delim");
  value_collection.def(pybind11::init<>());

  value_collection.def(pybind11::init([](const pybind11::dict& dict) -> ValueCollection {
                         ValueCollection collection;
                         update(collection, dict, false);
                         return collection;
                       }),
                       R"delim(
      Initialize the value collection from a python dictionary.

      >>> coll = ValueCollection({"a": 4, "b": 0.5})
      >>> len(coll)
      2
    )delim");

  value_collection.def("as_dict", &toDict, "Represent state as a Python dictionary");

  value_collection.def("update", &update, pybind11::arg("dict"), pybind11::arg("preserve_types") = true,
                       R"delim(
      Updates the collection with a dictionary

      :param dict: Dictionary to update the collection with
      :param preserve_types: Raise ``RuntimeError`` if an existing value would
        be overwritten by a different type

      >>> calculation_args = ValueCollection()
      >>> calculation_args["molecular_charge"] = 4
      >>> calculation_args["spin_multiplicity"] = 1
      >>> calculation_args.update({"molecular_charge": 0})
      >>> calculation_args["molecular_charge"]
      0
    )delim");

  value_collection.def(
      "get",
      [](const ValueCollection& coll, const std::string& key, pybind11::handle default_item) -> pybind11::handle {
        if (coll.valueExists(key)) {
          return pybind11::cast(GenericValueMeta::convert(coll.getValue(key)));
        }

        return default_item;
      },
      pybind11::arg("key"), pybind11::arg("default") = pybind11::none(),
      R"delim(
      Fetch a value from the collection by key. If the key is not in the
      collection, returns the default argument.

      :param key: Key to look up in the collection
      :param default: Value returned if the key is not in the map

      >>> coll = ValueCollection()
      >>> coll["a"] = 4
      >>> coll.get("a")
      4
      >>> coll.get("b") is None
      True
      >>> coll.get("b", 10)
      10
    )delim");

  value_collection.def("keys", &ValueCollection::getKeys, "List the keys in the collection");

  value_collection.def("__repr__",
                       [](const ValueCollection& coll) -> pybind11::str { return pybind11::str(toDict(coll)); });

  value_collection.def("__getitem__", [](const ValueCollection& coll, const std::string& key) -> GenericValueMeta::Variant {
    return GenericValueMeta::convert(coll.getValue(key));
  });

  value_collection.def("__setitem__", [](ValueCollection& coll, const std::string& key, const GenericValueMeta::Variant& v) {
    if (coll.valueExists(key)) {
      coll.modifyValue(key, GenericValueMeta::convert(v));
    }
    else {
      coll.addGenericValue(key, GenericValueMeta::convert(v));
    }
  });

  value_collection.def("__getitem__",
                       [](const ValueCollection& coll, const unsigned i) -> std::pair<std::string, GenericValueMeta::Variant> {
                         const auto key = coll.getKeys().at(i);
                         return std::make_pair(key, GenericValueMeta::convert(coll.getValue(key)));
                       });

  value_collection.def("__len__", &ValueCollection::count);

  value_collection.def("__contains__",
                       [](const ValueCollection& coll, const std::string& key) { return coll.valueExists(key); });
}
