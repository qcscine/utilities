/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Pybind.h"
#include "pybind11/operators.h"
#include "boost/bind.hpp"
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/GenericValueVariant.h>
#include <Utils/UniversalSettings/ParametrizedOptionValue.h>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/algorithm/iteration.hpp>
#include <boost/fusion/include/iteration.hpp>
#include <boost/fusion/include/std_tuple.hpp>

using namespace Scine::Utils;
using namespace UniversalSettings;

namespace {

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

    boost::optional<GenericValue> generic;
    try {
      generic = GenericValueMeta::convert(item.second.cast<GenericValueMeta::Variant>());
    }
    catch (pybind11::cast_error& /* e */) {
      const auto type = item.second.attr("__class__").attr("__name__").cast<std::string>();
      throw std::runtime_error("Could not store " + type + " at key '" + key + "': No matching GenericValue type");
    }
    if (coll.valueExists(key)) {
      if (preserveTypes && !GenericValueMeta::sameType(*generic, coll.getValue(key))) {
        throw std::runtime_error("GenericValue already exists as a different type!");
      }
      coll.modifyValue(key, std::move(*generic));
    }
    else {
      coll.addGenericValue(key, std::move(*generic));
    }
  }
}

struct CollectionIterWrapper {
  // Pybind doesn't need much by way of iterator requirements
  ValueCollection::Container::const_iterator iter;

  std::pair<std::string, GenericValueMeta::Variant> operator*() const {
    return std::make_pair(iter->first, GenericValueMeta::convert(iter->second));
  }

  CollectionIterWrapper& operator++() {
    ++iter;
    return *this;
  }

  bool operator==(const CollectionIterWrapper& other) const {
    return iter == other.iter;
  }
};

void init_value_collection(pybind11::module& m) {
  // "Forward-declare" ParametrizedOptionValue
  pybind11::class_<ParametrizedOptionValue> parametrizedOptionValue(m, "ParametrizedOptionValue");

  pybind11::class_<ValueCollection, std::shared_ptr<ValueCollection>> value_collection(m, "ValueCollection",
                                                                                       R"delim(
      Type-erased C++ map-like container with string keys that can store the
      following types of values: ``bool``, ``int``, ``float``, ``str``,
      ``ValueCollection`` (enables nesting!), ``List[int]``, ``List[str]``,
      and ``List[ValueCollection]``.

      Has members to imitate behavior of a python dictionary with string keys.

      >>> coll = ValueCollection()
      >>> coll
      scine_utilities.ValueCollection({})
      >>> coll["a"] = [1, 2, 3, 4]  # Add or lookup elements by key
      >>> coll
      scine_utilities.ValueCollection({'a': [1, 2, 3, 4]})
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

  value_collection.def(pybind11::init<const ValueCollection&>());

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
      [](const ValueCollection& coll, const std::string& key, pybind11::object default_item) -> pybind11::object {
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

  value_collection.def_static("from_dict", [](const pybind11::dict& dict) -> ValueCollection {
    ValueCollection collection;
    update(collection, dict, false);
    return collection;
  });

  value_collection.def("__repr__", [](pybind11::handle coll) -> std::string {
    return qualifiedName(coll) + "(" + pybind11::repr(toDict(coll.cast<ValueCollection>())).cast<std::string>() + ")";
  });

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

  value_collection.def("__len__", &ValueCollection::size);

  value_collection.def("__contains__",
                       [](const ValueCollection& coll, const std::string& key) { return coll.valueExists(key); });

  value_collection.def("__delitem__", [](ValueCollection& coll, const std::string& key) { coll.dropValue(key); });

  value_collection.def("items", [](const ValueCollection& coll) {
    return pybind11::make_iterator(CollectionIterWrapper{coll.begin()}, CollectionIterWrapper{coll.end()});
  });

  /* Pickling support would be better, but is very tricky seeing as
   * ValueCollections have nested class states, and using C++ copy constructors
   * might be preferable
   */
  value_collection.def("__deepcopy__", [](const ValueCollection& coll, const pybind11::dict & /* memo */) -> ValueCollection {
    return coll;
  });

  value_collection.def(pybind11::self == pybind11::self);
  value_collection.def(pybind11::self != pybind11::self);

  // And now for the parametrized option value methods
  parametrizedOptionValue.def(pybind11::init<std::string, ValueCollection>());
  parametrizedOptionValue.def_readwrite("selected_option", &ParametrizedOptionValue::selectedOption);
  parametrizedOptionValue.def_readwrite("option_settings", &ParametrizedOptionValue::optionSettings);

  parametrizedOptionValue.def(pybind11::self == pybind11::self);
  parametrizedOptionValue.def(pybind11::self != pybind11::self);
}
