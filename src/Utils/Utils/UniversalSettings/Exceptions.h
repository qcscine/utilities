/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_EXCEPTIONS_H
#define UNIVERSALSETTINGS_EXCEPTIONS_H
/* Internal Headers */
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <stdexcept>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

/*!
 * Base class for UniversalSettings exceptions.
 */
class Exception : public std::runtime_error {
 public:
  explicit Exception(const std::string& s) : std::runtime_error(s) {
  }
};

class InvalidDescriptorConversionException : public Exception {
 public:
  explicit InvalidDescriptorConversionException(const SettingDescriptor& d)
    : Exception("Error when trying to convert setting descriptor \"" + d.getPropertyDescription() + "\".") {
  }
};

class AlreadyExistingDescriptorException : public Exception {
 public:
  explicit AlreadyExistingDescriptorException(const std::string& name)
    : Exception("A GenericDescriptor with name \"" + name + "\" already exists in the DescriptorCollection.") {
  }
};

class InexistingDescriptorException : public Exception {
 public:
  explicit InexistingDescriptorException(const std::string& name)
    : Exception("No GenericDescriptor with name \"" + name + "\" exists in the DescriptorCollection.") {
  }
};

class InvalidValueConversionException : public Exception {
 public:
  InvalidValueConversionException() : Exception("Error when trying to convert a setting value.") {
  }
};

class AlreadyExistingValueException : public Exception {
 public:
  explicit AlreadyExistingValueException(const std::string& name)
    : Exception("A GenericValue with name \"" + name + "\" already exists in the ValueCollection.") {
  }
};

class InexistingValueException : public Exception {
 public:
  explicit InexistingValueException(const std::string& name)
    : Exception("No GenericValue with name \"" + name + "\" exists in the ValueCollection.") {
  }
};

class ValueHasDifferentTypeException : public Exception {
 public:
  explicit ValueHasDifferentTypeException(const std::string& propertyDescription)
    : Exception("The GenericValue \"" + propertyDescription + "\" has a different type than required.") {
  }
};

class EmptyOptionListException : public Exception {
 public:
  explicit EmptyOptionListException(const std::string& propertyDescription)
    : Exception("The OptionList property \"" + propertyDescription + "\" has no items.") {
  }
};

class OptionAlreadyExistsException : public Exception {
 public:
  OptionAlreadyExistsException(const std::string& option, const std::string& propertyDescription)
    : Exception("The OptionList property \"" + propertyDescription + "\" already contains an option called \"" + option + "\"") {
  }
};

class OptionDoesNotExistException : public Exception {
 public:
  OptionDoesNotExistException(const std::string& option, const std::string& propertyDescription)
    : Exception("The OptionList property \"" + propertyDescription + "\" contains no option called \"" + option + "\"") {
  }
};

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
#endif // UNIVERSALSETTINGS_EXCEPTIONS_H
