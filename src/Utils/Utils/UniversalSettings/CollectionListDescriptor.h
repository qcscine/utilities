/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_COLLECTIONLISTDESCRIPTOR_H
#define UNIVERSALSETTINGS_COLLECTIONLISTDESCRIPTOR_H
/* Internal Headers */
#include "Utils/UniversalSettings/DescriptorCollection.h"
#include "Utils/UniversalSettings/SettingDescriptor.h"
#include "Utils/UniversalSettings/ValueCollection.h"

namespace Scine {
namespace Utils {
namespace UniversalSettings {

/*!
 * @brief SettingDescriptor for one DescriptorCollection constraining a list of
 *   ValueCollections.
 */
class CollectionListDescriptor : public SettingDescriptor {
 public:
  CollectionListDescriptor() = default;
  /**
   * @brief Constructor
   *
   * @param propertyDescription String describing what the list of settings is for
   * @param base A list of SettingsDescriptors that set what 'fields' and
   *   default values exist for this setting
   */
  explicit CollectionListDescriptor(std::string propertyDescription, DescriptorCollection base);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;
  std::string explainInvalidValue(const GenericValue& v) const override;

  /**
   * @brief Checks whether a supplied list of string, GenericValues match the
   *   expected fields and valid values set by this instance
   *
   * @param v A vector<ValueCollection>, which in turn is just a list of string
   *   and GenericValue pairs, that map fields to values
   *
   * @return Whether the supplied object is valid with the configuration set
   */
  bool validValue(const GenericValue::CollectionList& v) const;

  /**
   * @brief Get the configuration of this setting
   *
   * @return The configuration object of this setting
   */
  const DescriptorCollection& getDescriptorCollection() const;

 private:
  DescriptorCollection base_;
};

inline CollectionListDescriptor::CollectionListDescriptor(std::string propertyDescription, DescriptorCollection base)
  : SettingDescriptor(std::move(propertyDescription)), base_(std::move(base)) {
}

inline std::unique_ptr<SettingDescriptor> CollectionListDescriptor::clone() const {
  return std::make_unique<CollectionListDescriptor>(*this);
}

inline GenericValue CollectionListDescriptor::getDefaultGenericValue() const {
  return GenericValue::fromCollectionList({});
}

inline bool CollectionListDescriptor::validValue(const GenericValue::CollectionList& v) const {
  return std::all_of(std::begin(v), std::end(v), [&](const auto& coll) { return base_.validValue(coll); });
}

inline bool CollectionListDescriptor::validValue(const GenericValue& v) const {
  if (!v.isCollectionList()) {
    return false;
  }
  return validValue(v.toCollectionList());
}

inline const DescriptorCollection& CollectionListDescriptor::getDescriptorCollection() const {
  return base_;
}

inline std::string CollectionListDescriptor::explainInvalidValue(const GenericValue& v) const {
  assert(!validValue(v));
  if (!v.isCollectionList()) {
    return "Generic value for collection list setting '" + getPropertyDescription() + "' is not a collection list!";
  }
  std::string explanation;
  for (const auto& coll : v.toCollectionList()) {
    explanation += base_.explainInvalidValue(coll);
  }
  return explanation;
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_COLLECTIONLISTDESCRIPTOR_H
