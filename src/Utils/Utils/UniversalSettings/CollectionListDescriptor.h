/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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
 * @brief SettingDescriptor for multiple SettingsDescriptors
 */
class CollectionListDescriptor : public SettingDescriptor {
 public:
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

  /**
   * @brief Checks whether a supplied list of string, GenericValues match the
   *   expected fields and valid values set by this instance
   *
   * @param v A vector<ValueCollection>, which in turn is just a list of string
   *   and GenericValue pairs, that map fields to values
   *
   * @return Whether the supplied object is valid with the configuration set
   */
  bool validValue(GenericValue::CollectionList v) const;

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

inline bool CollectionListDescriptor::validValue(GenericValue::CollectionList v) const {
  for (const auto& i : v) {
    if (!base_.validValue(i)) {
      return false;
    }
  }
  return true;
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

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_COLLECTIONLISTDESCRIPTOR_H
