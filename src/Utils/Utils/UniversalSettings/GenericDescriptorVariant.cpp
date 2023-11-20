/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/UniversalSettings/GenericDescriptorVariant.h"
#include "Utils/UniversalSettings/GenericValueVariant.h"
#include <Utils/UniversalSettings/CollectionListDescriptor.h>
#include <Utils/UniversalSettings/ParametrizedOptionListDescriptor.h>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/algorithm/iteration.hpp>
#include <boost/fusion/include/iteration.hpp>
#include <boost/fusion/include/std_tuple.hpp>

namespace Scine {
namespace Utils {
namespace UniversalSettings {
namespace {

template<typename VariantType, typename DescriptorType>
VariantType convertBase(DescriptorType& v) {
  using VariantOption = boost::optional<VariantType>;
  constexpr unsigned N = std::tuple_size<GenericDescriptorMeta::Types>::value;
  // clang-format off
  return boost::fusion::accumulate(
    GenericValueMeta::tuple_index_sequence<N>(),
    VariantOption(boost::none),
    [&](VariantOption maybeVariant, auto t) -> VariantOption {
      if (maybeVariant) {
        return maybeVariant;
      }

      constexpr unsigned i = decltype(t)::value;
      using T = std::tuple_element_t<i, GenericDescriptorMeta::Types>;
      if (v.template is<T>()) {
        return VariantType(v.template get<T>());
      }

      return maybeVariant;
    }
  ).value();
  // clang-format on
}

} // namespace

GenericDescriptorMeta::RefVariant GenericDescriptorMeta::convert(GenericDescriptor& v) {
  return convertBase<GenericDescriptorMeta::RefVariant>(v);
}

GenericDescriptorMeta::CRefVariant GenericDescriptorMeta::convert(const GenericDescriptor& v) {
  return convertBase<GenericDescriptorMeta::CRefVariant>(v);
}

} // namespace UniversalSettings
} // namespace Utils
} // namespace Scine
