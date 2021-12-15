/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_PARAMETRIZEDOPTIONVALUE_H
#define UNIVERSALSETTINGS_PARAMETRIZEDOPTIONVALUE_H
/* Internal Headers */
#include "Utils/UniversalSettings/ValueCollection.h"
/* External Headers */
#include <string>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

/*!
 * Value struct for a ParametrizedOptionList
 * NB: settingsKey is relevant for the generation of input files, as the option and the settings would be saved in
 * different nodes in YAML, for instance. The key for the selectedOption is not registered, it is taken to be
 * the one saved in SettingValueCollection.
 */
struct ParametrizedOptionValue {
  ParametrizedOptionValue(std::string option, ValueCollection settings);
  std::string selectedOption;
  ValueCollection optionSettings;

  bool operator==(const ParametrizedOptionValue& other) const {
    return (selectedOption == other.selectedOption && optionSettings == other.optionSettings);
  }

  bool operator!=(const ParametrizedOptionValue& other) const {
    return !(*this == other);
  }
};

inline ParametrizedOptionValue::ParametrizedOptionValue(std::string option, ValueCollection settings)
  : selectedOption(std::move(option)), optionSettings(std::move(settings)) {
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
#endif // UNIVERSALSETTINGS_PARAMETRIZEDOPTIONVALUE_H
