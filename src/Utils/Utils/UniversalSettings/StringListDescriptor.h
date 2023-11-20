/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_STRINGLISTDESCRIPTOR_H
#define UNIVERSALSETTINGS_STRINGLISTDESCRIPTOR_H
/* Internal Headers */
#include "GenericValue.h"
#include "SettingDescriptor.h"
/* External Headers */
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace Scine {
namespace Utils {
namespace UniversalSettings {
/*!
 * SettingDescriptor for a list of strings.
 */
class StringListDescriptor : public SettingDescriptor {
 public:
  using StringList = GenericValue::StringList;

  StringListDescriptor() = default;
  explicit StringListDescriptor(std::string propertyDescription);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;
  std::string explainInvalidValue(const GenericValue& v) const override;

  //! Returns the default string list value
  StringList getDefaultValue() const;
  //! Sets the default string list value
  void setDefaultValue(StringList def);

  //! Returns the default value for an individual string in the list
  std::string getDefaultItemValue() const;
  //! Get the default value for an individual string in the list
  void setDefaultItemValue(std::string def);

 private:
  StringList defaultValue_{};
  StringList::value_type defaultItemValue_{};
};

inline StringListDescriptor::StringListDescriptor(std::string propertyDescription)
  : SettingDescriptor(std::move(propertyDescription)) {
}

inline StringListDescriptor::StringList StringListDescriptor::getDefaultValue() const {
  return defaultValue_;
}

inline void StringListDescriptor::setDefaultValue(StringList def) {
  defaultValue_ = std::move(def);
}

inline std::unique_ptr<SettingDescriptor> StringListDescriptor::clone() const {
  return std::make_unique<StringListDescriptor>(*this);
}

inline GenericValue StringListDescriptor::getDefaultGenericValue() const {
  return GenericValue::fromStringList(getDefaultValue());
}

inline bool StringListDescriptor::validValue(const GenericValue& v) const {
  return v.isStringList();
}

inline std::string StringListDescriptor::getDefaultItemValue() const {
  return defaultItemValue_;
}

inline void StringListDescriptor::setDefaultItemValue(std::string def) {
  defaultItemValue_ = std::move(def);
}

inline std::string StringListDescriptor::explainInvalidValue(const GenericValue& /* v */) const {
  return "Generic value for string setting '" + getPropertyDescription() + "' is not a string list!";
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_STRINGLISTDESCRIPTOR_H
