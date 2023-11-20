/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_DIRECTORYDESCRIPTOR_H
#define UNIVERSALSETTINGS_DIRECTORYDESCRIPTOR_H
/* Internal Headers */
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <string>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

/*!
 * SettingDescriptor for a directory path.
 */
class DirectoryDescriptor : public SettingDescriptor {
 public:
  DirectoryDescriptor() = default;
  /**
   * @brief Constructor
   *
   * @param propertyDescription The string describing what the directory path is for
   */
  explicit DirectoryDescriptor(std::string propertyDescription);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;
  std::string explainInvalidValue(const GenericValue& v) const override;

  //! Gets the default value of the directory path
  const std::string& getDefaultValue() const;

  //! Sets the default value of the directory path
  void setDefaultValue(std::string def);

 private:
  std::string defaultValue_{};
};

inline DirectoryDescriptor::DirectoryDescriptor(std::string propertyDescription)
  : SettingDescriptor(std::move(propertyDescription)) {
}

inline const std::string& DirectoryDescriptor::getDefaultValue() const {
  return defaultValue_;
}

inline void DirectoryDescriptor::setDefaultValue(std::string def) {
  defaultValue_ = std::move(def);
}

inline std::unique_ptr<SettingDescriptor> DirectoryDescriptor::clone() const {
  return std::make_unique<DirectoryDescriptor>(*this);
}

inline GenericValue DirectoryDescriptor::getDefaultGenericValue() const {
  return GenericValue::fromString(getDefaultValue());
}

inline bool DirectoryDescriptor::validValue(const GenericValue& v) const {
  return v.isString();
}

inline std::string DirectoryDescriptor::explainInvalidValue(const GenericValue& /* v */) const {
  return "Generic value for string setting '" + getPropertyDescription() + "' is not a string!";
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_DIRECTORYDESCRIPTOR_H
