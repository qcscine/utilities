/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_SETTINGDESCRIPTOR_H
#define UNIVERSALSETTINGS_SETTINGDESCRIPTOR_H
/* External Headers */
#include <memory>
#include <string>

namespace Scine {
namespace Utils {
namespace UniversalSettings {
/* Forward Declarations */
class GenericValue;

/*!
 * @brief Class to restrict a particular setting's type and values.
 *
 * Contains a string explaining what the setting is for, possibly a
 * default value and a function that checks if the contained value is valid for
 * the setting described.
 *
 * This is an abstract base class for type-specific descriptor implementations.
 * Other classes derive from it not because of the need for polymorphism, but
 * to collect common functionality (i.e. property description).
 */
class SettingDescriptor {
 public:
  virtual ~SettingDescriptor() = default;

  /**
   * @brief Create a heap-allocated base pointer copy of this instance
   *
   * @return A heap-allocated copy of this instance
   */
  virtual std::unique_ptr<SettingDescriptor> clone() const = 0;

  /**
   * @brief Changes the string explanation of what the setting is for
   *
   * @param propertyDescription The string explanation of what the setting is for
   */
  void setPropertyDescription(std::string propertyDescription);

  /**
   * @brief Fetches the string explanation of what the setting is for
   *
   * @return The string explanation of what the setting is for
   */
  const std::string& getPropertyDescription() const;

  /**
   * @brief Creates a type-erased representation of the setting's default value
   *
   * @return A type-erased representation of the setting's default value
   */
  virtual GenericValue getDefaultGenericValue() const = 0;

  /**
   * @brief Checks if a particular type-erased representation of a setting is a
   *   valid value for this kind of setting.
   *
   * @param v A type-erased setting representation
   *
   * @return Whether the type-erased representation is a valid value for this
   *   kind of setting
   */
  virtual bool validValue(const GenericValue& v) const = 0;

 protected:
  /*! Allow instantiation only through derived classes. */
  explicit SettingDescriptor(std::string propertyDescription);

 private:
  /**
   * @brief String that explains what the setting is for
   */
  std::string propertyDescription_;
};

inline SettingDescriptor::SettingDescriptor(std::string propertyDescription)
  : propertyDescription_(std::move(propertyDescription)) {
}

inline void SettingDescriptor::setPropertyDescription(std::string propertyDescription) {
  propertyDescription_ = std::move(propertyDescription);
}

inline const std::string& SettingDescriptor::getPropertyDescription() const {
  return propertyDescription_;
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
#endif // UNIVERSALSETTINGS_SETTINGDESCRIPTOR_H
