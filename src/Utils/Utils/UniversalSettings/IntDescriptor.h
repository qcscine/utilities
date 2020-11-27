/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_INTDESCRIPTOR_H
#define UNIVERSALSETTINGS_INTDESCRIPTOR_H
/* Internal Headers */
#include "Utils/UniversalSettings/GenericValue.h"
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <limits>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

/*!
 * SettingDescriptor for an integer value.
 */
class IntDescriptor : public SettingDescriptor {
 public:
  /**
   * @brief Constructor
   *
   * @param propertyDescription A string describing what the integer is for
   */
  explicit IntDescriptor(std::string propertyDescription);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;

  //! Checks whether the supplied value is within the configured bounds
  bool validValue(int v) const;

  //! Returns the lower bound on valid values
  int getMinimum() const;
  //! Returns the upper bound on valid values
  int getMaximum() const;
  //! Returns the default value
  int getDefaultValue() const;

  //! Sets the lower bound on valid values
  void setMinimum(int min);
  //! Sets the upper bound on valid values
  void setMaximum(int max);
  //! Sets the default value
  void setDefaultValue(int def);

 private:
  int minimum_;
  int maximum_;
  int defaultValue_;
};

inline IntDescriptor::IntDescriptor(std::string propertyDescription)
  : SettingDescriptor(std::move(propertyDescription)) {
  minimum_ = std::numeric_limits<int>::lowest();
  maximum_ = std::numeric_limits<int>::max();
  defaultValue_ = 0;
}

inline int IntDescriptor::getMinimum() const {
  return minimum_;
}

inline int IntDescriptor::getMaximum() const {
  return maximum_;
}

inline int IntDescriptor::getDefaultValue() const {
  return defaultValue_;
}

inline void IntDescriptor::setMinimum(int min) {
  minimum_ = min;
}

inline void IntDescriptor::setMaximum(int max) {
  maximum_ = max;
}

inline void IntDescriptor::setDefaultValue(int def) {
  defaultValue_ = def;
}

inline std::unique_ptr<SettingDescriptor> IntDescriptor::clone() const {
  return std::make_unique<IntDescriptor>(*this);
}

inline GenericValue IntDescriptor::getDefaultGenericValue() const {
  return GenericValue::fromInt(getDefaultValue());
}

inline bool IntDescriptor::validValue(int v) const {
  return getMinimum() <= v && v <= getMaximum();
}

inline bool IntDescriptor::validValue(const GenericValue& v) const {
  if (!v.isInt()) {
    return false;
  };
  return validValue(v.toInt());
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_INTDESCRIPTOR_H
