/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INCLUDE_UTILS_UNIVERSAL_SETTINGS_GENERIC_DESCRIPTOR_VARIANT_H
#define INCLUDE_UTILS_UNIVERSAL_SETTINGS_GENERIC_DESCRIPTOR_VARIANT_H

#include "Utils/UniversalSettings/GenericDescriptor.h"
#include "boost/variant.hpp"
#include <functional>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

struct GenericDescriptorMeta {
  using Types = std::tuple<BoolDescriptor, IntDescriptor, DoubleDescriptor, StringDescriptor, FileDescriptor,
                           DirectoryDescriptor, OptionListDescriptor, DescriptorCollection, ParametrizedOptionListDescriptor,
                           IntListDescriptor, DoubleListDescriptor, StringListDescriptor, CollectionListDescriptor>;

  template<typename Tuple>
  struct RewriteAsVariantOfRefs {
    template<typename std::size_t... Inds>
    static constexpr auto impl(std::index_sequence<Inds...> inds) -> boost::variant<std::tuple_element_t<Inds, Tuple>&...>;

    using Type = decltype(impl(std::make_index_sequence<std::tuple_size<Tuple>::value>()));
  };

  template<typename Tuple>
  struct RewriteAsVariantOfCRefs {
    template<typename std::size_t... Inds>
    static constexpr auto impl(std::index_sequence<Inds...> inds)
        -> boost::variant<const std::tuple_element_t<Inds, Tuple>&...>;

    using Type = decltype(impl(std::make_index_sequence<std::tuple_size<Tuple>::value>()));
  };

  using RefVariant = RewriteAsVariantOfRefs<Types>::Type;
  using CRefVariant = RewriteAsVariantOfCRefs<Types>::Type;

  //! Converts a generic descriptor into a variant of references
  static RefVariant convert(GenericDescriptor& v);
  //! Converts a generic descriptor into a variant of const references
  static CRefVariant convert(const GenericDescriptor& v);
};

} // namespace UniversalSettings
} // namespace Utils
} // namespace Scine

#endif
