/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_INTLISTLISTDESCRIPTOR_H
#define UTILSOS_INTLISTLISTDESCRIPTOR_H
/* Internal Headers */
#include "Utils/UniversalSettings/Exceptions.h"
#include "Utils/UniversalSettings/GenericValue.h"
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <algorithm>
#include <limits>
#include <vector>

namespace Scine {
namespace Utils {
namespace UniversalSettings {
/*!
 * SettingDescriptor for a list of integer lists.
 */
class IntListListDescriptor : public SettingDescriptor {
 public:
  using IntListList = GenericValue::IntListList;

  IntListListDescriptor() = default;
  /**
   * @brief Constructor
   *
   * @param propertyDescription A string describing what the list of integer lists is for
   */
  explicit IntListListDescriptor(std::string propertyDescription);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;
  std::string explainInvalidValue(const GenericValue& v) const override;

  // Checks that all contained integers are valid
  bool validValue(const IntListList& v) const;

  //! Get the default list of integers
  IntListList getDefaultValue() const;
  //! Set the default list of integers
  void setDefaultValue(IntListList def);

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

  IntListList defaultValue_{};
};

inline IntListListDescriptor::IntListListDescriptor(std::string propertyDescription)
  : SettingDescriptor(std::move(propertyDescription)) {
  itemMinimum_ = std::numeric_limits<int>::lowest();
  itemMaximum_ = std::numeric_limits<int>::max();
  defaultItemValue_ = 0;
}

inline IntListListDescriptor::IntListList IntListListDescriptor::getDefaultValue() const {
  return defaultValue_;
}

inline void IntListListDescriptor::setDefaultValue(IntListListDescriptor::IntListList def) {
  defaultValue_ = std::move(def);
}

inline int IntListListDescriptor::getItemMinimum() const {
  return itemMinimum_;
}

inline int IntListListDescriptor::getItemMaximum() const {
  return itemMaximum_;
}

inline int IntListListDescriptor::getDefaultItemValue() const {
  return defaultItemValue_;
}

inline void IntListListDescriptor::setItemMinimum(int min) {
  itemMinimum_ = min;
}

inline void IntListListDescriptor::setItemMaximum(int max) {
  itemMaximum_ = max;
}

inline void IntListListDescriptor::setDefaultItemValue(int def) {
  defaultItemValue_ = def;
}

inline std::unique_ptr<SettingDescriptor> IntListListDescriptor::clone() const {
  return std::make_unique<IntListListDescriptor>(*this);
}

inline GenericValue IntListListDescriptor::getDefaultGenericValue() const {
  return GenericValue::fromIntListList(getDefaultValue());
}

inline bool IntListListDescriptor::validValue(const IntListList& v) const {
  const int minimum = getItemMinimum();
  const int maximum = getItemMaximum();
  for (const auto& list : v) {
    if (!std::all_of(std::begin(list), std::end(list), [=](int i) { return minimum <= i && i <= maximum; })) {
      return false;
    }
  }
  return true;
}

inline bool IntListListDescriptor::validValue(const GenericValue& v) const {
  if (!v.isIntListList()) {
    return false;
  }
  return validValue(v.toIntListList());
}

inline std::string IntListListDescriptor::explainInvalidValue(const GenericValue& v) const {
  assert(!validValue(v));
  if (!v.isIntListList()) {
    return "Generic value for integer list setting '" + getPropertyDescription() + "' is not a list of integer lists!";
  }
  std::string explanation = "A value in the list of integer lists descriptor '" + getPropertyDescription() + "' is ";
  explanation += "out of bounds [" + std::to_string(itemMinimum_) + ", " + std::to_string(itemMaximum_) + "].";
  return explanation;
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
#endif // UTILSOS_INTLISTLISTDESCRIPTOR_H
