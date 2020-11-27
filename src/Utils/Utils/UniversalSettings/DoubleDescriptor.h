/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_DOUBLEDESCRIPTOR_H
#define UNIVERSALSETTINGS_DOUBLEDESCRIPTOR_H
/* Internal Headers */
#include "Utils/UniversalSettings/GenericValue.h"
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <limits>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

/*!
 * @brief SettingDescriptor for a floating-point value.
 *
 * Can set boundaries on valid values and a default value for this kind of
 * setting.
 */
class DoubleDescriptor : public SettingDescriptor {
 public:
  /**
   * @brief Constructor
   *
   * @param propertyDescription A string describing what the double is for
   */
  explicit DoubleDescriptor(std::string propertyDescription);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;

  //! Checks whether the supplied value is within the configured bounds
  bool validValue(double v) const;

  //! Returns the lower bound on valid values
  double getMinimum() const;
  //! Returns the upper bound on valid values
  double getMaximum() const;
  //! Returns the default value
  double getDefaultValue() const;

  //! Sets the lower bound on valid values
  void setMinimum(double min);
  //! Sets the upper bound on valid values
  void setMaximum(double max);
  //! Sets the default value
  void setDefaultValue(double def);

 private:
  double minimum_;
  double maximum_;
  double defaultValue_;
};

inline DoubleDescriptor::DoubleDescriptor(std::string propertyDescription)
  : SettingDescriptor(std::move(propertyDescription)) {
  minimum_ = std::numeric_limits<double>::lowest();
  maximum_ = std::numeric_limits<double>::max();
  defaultValue_ = 0.0;
}

inline double DoubleDescriptor::getMinimum() const {
  return minimum_;
}

inline double DoubleDescriptor::getMaximum() const {
  return maximum_;
}

inline double DoubleDescriptor::getDefaultValue() const {
  return defaultValue_;
}

inline void DoubleDescriptor::setMinimum(double min) {
  minimum_ = min;
}

inline void DoubleDescriptor::setMaximum(double max) {
  maximum_ = max;
}

inline void DoubleDescriptor::setDefaultValue(double def) {
  defaultValue_ = def;
}

inline std::unique_ptr<SettingDescriptor> DoubleDescriptor::clone() const {
  return std::make_unique<DoubleDescriptor>(*this);
}

inline GenericValue DoubleDescriptor::getDefaultGenericValue() const {
  return GenericValue::fromDouble(getDefaultValue());
}

inline bool DoubleDescriptor::validValue(const GenericValue& v) const {
  if (!v.isDouble()) {
    return false;
  }
  return validValue(v.toDouble());
}

inline bool DoubleDescriptor::validValue(double v) const {
  return getMinimum() <= v && v <= getMaximum();
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_DOUBLEDESCRIPTOR_H
