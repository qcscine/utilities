/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_DOUBLELISTDESCRIPTOR_H
#define UNIVERSALSETTINGS_DOUBLELISTDESCRIPTOR_H
/* Doubleernal Headers */
#include "Utils/UniversalSettings/Exceptions.h"
#include "Utils/UniversalSettings/GenericValue.h"
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <algorithm>
#include <limits>
#include <sstream>
#include <vector>

namespace Scine {
namespace Utils {
namespace UniversalSettings {
/*!
 * SettingDescriptor for a list of double values.
 */
class DoubleListDescriptor : public SettingDescriptor {
 public:
  using DoubleList = GenericValue::DoubleList;

  DoubleListDescriptor() = default;
  /**
   * @brief Constructor
   *
   * @param propertyDescription A string describing what the list of doubles is for
   */
  explicit DoubleListDescriptor(std::string propertyDescription);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;
  std::string explainInvalidValue(const GenericValue& v) const override;

  // Checks that all contained doubles are valid
  bool validValue(const DoubleList& v) const;

  //! Get the default list of doubles
  DoubleList getDefaultValue() const;
  //! Set the default list of doubles
  void setDefaultValue(DoubleList def);

  //! Returns the default value for an item in the doubles list
  double getDefaultItemValue() const;
  //! Set the default value for each item in the doubles list
  void setDefaultItemValue(double def);

  //! Returns the lower bound for each value in the doubles list
  double getItemMinimum() const;
  //! Set the lower bound for each value in the doubles list
  void setItemMinimum(double min);

  //! Returns the upper bound for each value in the doubles list
  double getItemMaximum() const;
  //! Set the upper bound for each value in the doubles list
  void setItemMaximum(double max);

 private:
  double itemMinimum_;
  double itemMaximum_;
  double defaultItemValue_;

  DoubleList defaultValue_{};
};

inline DoubleListDescriptor::DoubleListDescriptor(std::string propertyDescription)
  : SettingDescriptor(std::move(propertyDescription)) {
  itemMinimum_ = std::numeric_limits<double>::lowest();
  itemMaximum_ = std::numeric_limits<double>::max();
  defaultItemValue_ = 0;
}

inline DoubleListDescriptor::DoubleList DoubleListDescriptor::getDefaultValue() const {
  return defaultValue_;
}

inline void DoubleListDescriptor::setDefaultValue(DoubleListDescriptor::DoubleList def) {
  defaultValue_ = std::move(def);
}

inline double DoubleListDescriptor::getItemMinimum() const {
  return itemMinimum_;
}

inline double DoubleListDescriptor::getItemMaximum() const {
  return itemMaximum_;
}

inline double DoubleListDescriptor::getDefaultItemValue() const {
  return defaultItemValue_;
}

inline void DoubleListDescriptor::setItemMinimum(double min) {
  itemMinimum_ = min;
}

inline void DoubleListDescriptor::setItemMaximum(double max) {
  itemMaximum_ = max;
}

inline void DoubleListDescriptor::setDefaultItemValue(double def) {
  defaultItemValue_ = def;
}

inline std::unique_ptr<SettingDescriptor> DoubleListDescriptor::clone() const {
  return std::make_unique<DoubleListDescriptor>(*this);
}

inline GenericValue DoubleListDescriptor::getDefaultGenericValue() const {
  return GenericValue::fromDoubleList(getDefaultValue());
}

inline bool DoubleListDescriptor::validValue(const DoubleList& v) const {
  const double minimum = getItemMinimum();
  const double maximum = getItemMaximum();

  return std::all_of(std::begin(v), std::end(v), [=](double i) { return minimum <= i && i <= maximum; });
}

inline bool DoubleListDescriptor::validValue(const GenericValue& v) const {
  if (!v.isDoubleList()) {
    return false;
  }
  return validValue(v.toDoubleList());
}

inline std::string DoubleListDescriptor::explainInvalidValue(const GenericValue& v) const {
  assert(!validValue(v));
  if (!v.isDoubleList()) {
    return "Generic value for double list setting '" + getPropertyDescription() + "' is not a double list!";
  }
  std::ostringstream explanation;
  explanation << "A value in the double list descriptor '" + getPropertyDescription() + "' is out of bounds ["
              << itemMinimum_ << "," << itemMaximum_ << "].";
  return explanation.str();
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
#endif // UNIVERSALSETTINGS_DOUBLELISTDESCRIPTOR
