/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_STRINGDESCRIPTOR_H
#define UNIVERSALSETTINGS_STRINGDESCRIPTOR_H
/* Internal Headers */
#include "Utils/UniversalSettings/GenericValue.h"
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <memory>
#include <string>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

/*!
 * SettingDescriptor for a string value.
 */
class StringDescriptor : public SettingDescriptor {
 public:
  StringDescriptor() = default;
  /**
   * @brief Constructor
   *
   * @param propertyDescription A string describing what the string setting is for
   */
  explicit StringDescriptor(std::string propertyDescription);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;
  std::string explainInvalidValue(const GenericValue& v) const override;

  //! Returns the default value
  const std::string& getDefaultValue() const;

  //! Sets the default string value
  void setDefaultValue(std::string def);

 private:
  std::string defaultValue_{};
};

inline StringDescriptor::StringDescriptor(std::string propertyDescription)
  : SettingDescriptor(std::move(propertyDescription)) {
}

inline const std::string& StringDescriptor::getDefaultValue() const {
  return defaultValue_;
}

inline void StringDescriptor::setDefaultValue(std::string def) {
  defaultValue_ = std::move(def);
}

inline std::unique_ptr<SettingDescriptor> StringDescriptor::clone() const {
  return std::make_unique<StringDescriptor>(*this);
}

inline GenericValue StringDescriptor::getDefaultGenericValue() const {
  return GenericValue::fromString(getDefaultValue());
}

inline bool StringDescriptor::validValue(const GenericValue& v) const {
  return v.isString();
}

inline std::string StringDescriptor::explainInvalidValue(const GenericValue& /* v */) const {
  return "Generic value for string setting '" + getPropertyDescription() + "' is not a string!";
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_STRINGDESCRIPTOR_H
