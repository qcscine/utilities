/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/UniversalSettings/GenericValueVariant.h"
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/algorithm/iteration.hpp>
#include <boost/fusion/include/iteration.hpp>
#include <boost/fusion/include/std_tuple.hpp>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

GenericValueMeta::Variant GenericValueMeta::convert(const GenericValue& v) {
  using VariantOption = boost::optional<Variant>;
  // clang-format off
  auto maybeVariant = boost::fusion::accumulate(
    GenericValueMeta::zip(
      GenericValueMeta::type_checkers(),
      GenericValueMeta::getters()
    ),
    VariantOption(boost::none),
    [&](VariantOption maybeVariant, const auto& t) -> VariantOption {
      if (maybeVariant) {
        return maybeVariant;
      }

      auto&& check_fn = std::get<0>(t);
      auto&& getter_fn = std::get<1>(t);
      if (check_fn(v)) {
        return Variant(getter_fn(v));
      }

      return maybeVariant;
    }
  );
  // clang-format on
  return maybeVariant.value();
}

GenericValue GenericValueMeta::convert(const Variant& v) {
  constexpr unsigned N = std::tuple_size<GenericValueMeta::Types>::value;
  using ValueOption = boost::optional<GenericValue>;
  // clang-format off
  auto maybeGenericValue = boost::fusion::accumulate(
    GenericValueMeta::zip(
      GenericValueMeta::tuple_index_sequence<N>(),
      GenericValueMeta::factories()
    ),
    ValueOption(boost::none),
    [&](ValueOption maybeGenericValue, auto t) -> ValueOption {
      if(maybeGenericValue) {
        return maybeGenericValue;
      }

      const auto type_index_constant = std::get<0>(t);
      auto&& factory_fn = std::get<1>(t);

      constexpr std::size_t variantTypeIndex = decltype(type_index_constant)::value;
      if (v.which() == variantTypeIndex) {
        using Subtype = std::tuple_element_t<variantTypeIndex, Types>;
        return factory_fn(boost::get<Subtype>(v));
      }

      return maybeGenericValue;
    }
  );
  // clang-format on
  return maybeGenericValue.value();
}

bool GenericValueMeta::sameType(const GenericValue& a, const GenericValue& b) {
  return boost::fusion::accumulate(GenericValueMeta::type_checkers(), false, [&](bool sameType, auto&& typecheckFn) -> bool {
    return sameType || (typecheckFn(a) && typecheckFn(b));
  });
}

} // namespace UniversalSettings
} // namespace Utils
} // namespace Scine
