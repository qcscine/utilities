/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Pybind.h>

using namespace Scine::Utils;

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
void recursivePropertyAdditionToResults(pybind11::class_<Results>& /* results */) {
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
