/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_INTLISTDESCRIPTOR_H
#define UNIVERSALSETTINGS_INTLISTDESCRIPTOR_H
/* Internal Headers */
#include "Utils/UniversalSettings/Exceptions.h"
#include "Utils/UniversalSettings/GenericValue.h"
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <limits>
#include <vector>

namespace Scine {
namespace Utils {
namespace UniversalSettings {
/*!
 * SettingDescriptor for a list of integer values.
 */
class IntListDescriptor : public SettingDescriptor {
 public:
  using IntList = GenericValue::IntList;

  /**
   * @brief Constructor
   *
   * @param propertyDescription A string describing what the list of integers is for
   */
  explicit IntListDescriptor(std::string propertyDescription);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;

  // Checks that all contained integers are valid
  bool validValue(IntList v) const;

  //! Get the default list of integers
  IntList getDefaultValue() const;
  //! Set the default list of integers
  void setDefaultValue(IntList def);

  //! Returns the default value for an item in the integer list
  int getDefaultItemValue() const;
  //! Set the default value for each item in the integer list
  void setDefaultItemValue(int def);

  //! Returns the lower bound for each value in the integer list
  int getItemMinimum() const;
  //! Set the lower bound for each value in the integer list
  void setItemMinimum(int min);

  //! Returns the upper bound for each value in the integer list
  int getItemMaximum() const;
  //! Set the upper bound for each value in the integer list
  void setItemMaximum(int max);

 private:
  int itemMinimum_;
  int itemMaximum_;
  int defaultItemValue_;

  IntList defaultValue_{};
};

inline IntListDescriptor::IntListDescriptor(std::string propertyDescription)
  : SettingDescriptor(std::move(propertyDescription)) {
  itemMinimum_ = std::numeric_limits<int>::lowest();
  itemMaximum_ = std::numeric_limits<int>::max();
  defaultItemValue_ = 0;
}

inline IntListDescriptor::IntList IntListDescriptor::getDefaultValue() const {
  return defaultValue_;
}

inline void IntListDescriptor::setDefaultValue(IntListDescriptor::IntList def) {
  defaultValue_ = std::move(def);
}

inline int IntListDescriptor::getItemMinimum() const {
  return itemMinimum_;
}

inline int IntListDescriptor::getItemMaximum() const {
  return itemMaximum_;
}

inline int IntListDescriptor::getDefaultItemValue() const {
  return defaultItemValue_;
}

inline void IntListDescriptor::setItemMinimum(int min) {
  itemMinimum_ = min;
}

inline void IntListDescriptor::setItemMaximum(int max) {
  itemMaximum_ = max;
}

inline void IntListDescriptor::setDefaultItemValue(int def) {
  defaultItemValue_ = def;
}

inline std::unique_ptr<SettingDescriptor> IntListDescriptor::clone() const {
  return std::make_unique<IntListDescriptor>(*this);
}

inline GenericValue IntListDescriptor::getDefaultGenericValue() const {
  return GenericValue::fromIntList(getDefaultValue());
}

inline bool IntListDescriptor::validValue(IntList v) const {
  for (const auto& i : v) {
    if (i < getItemMinimum() || getItemMaximum() < i) {
      return false;
    }
  }
  return true;
}

inline bool IntListDescriptor::validValue(const GenericValue& v) const {
  if (!v.isIntList()) {
    return false;
  }
  return validValue(v.toIntList());
}
} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
#endif // UNIVERSALSETTINGS_INTLISTDESCRIPTOR
