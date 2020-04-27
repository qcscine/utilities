/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/CalculatorBasics/Results.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/optional.hpp>

using namespace Scine::Utils;

// boost optional type caster
namespace pybind11 {
namespace detail {
template<typename T>
struct type_caster<boost::optional<T>> : optional_caster<boost::optional<T>> {};
} // namespace detail
} // namespace pybind11

namespace detail {

template<Property property>
void add_property(pybind11::class_<Results>& pyclass_results) {
  using ReturnType = std::tuple_element_t<getPropertyIndex(property), PropertyTypeTuple>;
  pyclass_results.def_property(
      propertyTypeName(property),
      [=](Results& results) -> boost::optional<ReturnType> {
        if (!std::mem_fn(&Results::has<property>)(results)) {
          return boost::none;
        }

        return std::mem_fn(&Results::get<property>)(results);
      },
      [=](Results& results, const boost::optional<ReturnType>& valueOptional) -> void {
        if (valueOptional) {
          std::mem_fn (&Results::set<property>)(results, *valueOptional);
        }
      });
}

template<typename T = void>
void recursivePropertyAdditionToResults(pybind11::class_<Results>& results) {
}

template<std::size_t index, std::size_t... Inds>
void recursivePropertyAdditionToResults(pybind11::class_<Results>& results) {
  add_property<static_cast<Property>(1 << index)>(results);

  recursivePropertyAdditionToResults<Inds...>(results);
}

template<std::size_t... Inds>
void unpackVariadicIndices(pybind11::class_<Results>& results, std::index_sequence<Inds...> /* inds */) {
  recursivePropertyAdditionToResults<Inds...>(results);
}

} // namespace detail

void init_results(pybind11::module& m) {
  pybind11::class_<Results> results(m, "Results");
  results.def(pybind11::init<>());

  ::detail::unpackVariadicIndices(results, std::make_index_sequence<std::tuple_size<PropertyTypeTuple>::value>{});
}
