/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_BOOLDESCRIPTOR_H
#define UNIVERSALSETTINGS_BOOLDESCRIPTOR_H
/* Internal Headers */
#include "Utils/UniversalSettings/GenericValue.h"
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <memory>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

/*! @brief Setting descriptor-derived boolean value
 *
 * SettingDescriptor for a bool value. This is a class that contains a string
 * that explains what the setting is for and also a default value.
 */
class BoolDescriptor : public SettingDescriptor {
 public:
  //! Construct with string describing what the value is for
  explicit BoolDescriptor(std::string propertyDescription);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;

  /**
   * @brief Get the default value
   *
   * @return Default value
   */
  bool getDefaultValue() const;

  /**
   * @brief Set the default value
   *
   * @param def The value to be set as default
   */
  void setDefaultValue(bool def);

 private:
  bool defaultValue_ = false;
};

inline BoolDescriptor::BoolDescriptor(std::string propertyDescription)
  : SettingDescriptor(std::move(propertyDescription)) {
}

inline std::unique_ptr<SettingDescriptor> BoolDescriptor::clone() const {
  return std::make_unique<BoolDescriptor>(*this);
}

inline GenericValue BoolDescriptor::getDefaultGenericValue() const {
  return GenericValue::fromBool(defaultValue_);
}

inline bool BoolDescriptor::getDefaultValue() const {
  return defaultValue_;
}

inline void BoolDescriptor::setDefaultValue(bool def) {
  defaultValue_ = def;
}

inline bool BoolDescriptor::validValue(const GenericValue& v) const {
  return v.isBool();
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_BOOLDESCRIPTOR_H
