/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_EXCEPTIONS_H
#define UNIVERSALSETTINGS_EXCEPTIONS_H
/* Internal Headers */
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <map>
#include <stdexcept>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

/*!
 * Base class for UniversalSettings exceptions.
 */
struct Exception : public std::runtime_error {
  explicit Exception(const std::string& s) : std::runtime_error(s) {
  }
};

struct InvalidDescriptorConversionException : public Exception {
  explicit InvalidDescriptorConversionException(const SettingDescriptor& d)
    : Exception("Error when trying to convert setting descriptor \"" + d.getPropertyDescription() + "\".") {
  }
};

struct AlreadyExistingDescriptorException : public Exception {
  explicit AlreadyExistingDescriptorException(const std::string& name)
    : Exception("A GenericDescriptor with name \"" + name + "\" already exists in the DescriptorCollection.") {
  }
};

struct InexistingDescriptorException : public Exception {
  explicit InexistingDescriptorException(const std::string& name)
    : Exception("No GenericDescriptor with name \"" + name + "\" exists in the DescriptorCollection.") {
  }
};

struct InvalidValueConversionException : public Exception {
  InvalidValueConversionException() : Exception("Error when trying to convert a setting value.") {
  }
};

struct InvalidSettingsException : public Exception {
  InvalidSettingsException(const std::string& explanation) : Exception(explanation) {
  }
};

struct AlreadyExistingValueException : public Exception {
  explicit AlreadyExistingValueException(const std::string& name)
    : Exception("A GenericValue with name \"" + name + "\" already exists in the ValueCollection.") {
  }
};

struct InexistingValueException : public Exception {
  explicit InexistingValueException(const std::string& name)
    : Exception("No GenericValue with name \"" + name + "\" exists in the ValueCollection.") {
  }
};

struct ValueHasDifferentTypeException : public Exception {
  explicit ValueHasDifferentTypeException(const std::string& propertyDescription)
    : Exception("The GenericValue \"" + propertyDescription + "\" has a different type than required.") {
  }
};

struct EmptyOptionListException : public Exception {
  explicit EmptyOptionListException(const std::string& propertyDescription)
    : Exception("The OptionList property \"" + propertyDescription + "\" has no items.") {
  }
};

struct OptionAlreadyExistsException : public Exception {
  OptionAlreadyExistsException(const std::string& option, const std::string& propertyDescription)
    : Exception("The OptionList property \"" + propertyDescription + "\" already contains an option called \"" + option + "\"") {
  }
};

struct OptionDoesNotExistException : public Exception {
  OptionDoesNotExistException(const std::string& option, const std::string& propertyDescription)
    : Exception("The OptionList property \"" + propertyDescription + "\" contains no option called \"" + option + "\"") {
  }
};

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
#endif // UNIVERSALSETTINGS_EXCEPTIONS_H
