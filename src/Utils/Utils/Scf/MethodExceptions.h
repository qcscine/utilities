/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_METHODEXCEPTIONS_H
#define UTILS_METHODEXCEPTIONS_H

#include "Utils/Geometry/ElementInfo.h"
#include <stdexcept>
#include <string>

namespace Scine {
namespace Utils {

namespace Methods {

class InitializationException : public std::runtime_error {
 public:
  explicit InitializationException(const std::string& what_arg) : runtime_error(what_arg) {
  }
};

class ParameterFileCannotBeOpenedException : public InitializationException {
 public:
  explicit ParameterFileCannotBeOpenedException(const std::string& filePath)
    : InitializationException("The following parameter file cannot be opened: " + filePath) {
  }
};

class ParameterFileIsInvalidException : public InitializationException {
 public:
  explicit ParameterFileIsInvalidException(const std::string& filePath)
    : InitializationException("The following parameter file is invalid: " + filePath) {
  }
};

class ParametersDoNotExistForElementException : public InitializationException {
 public:
  explicit ParametersDoNotExistForElementException(Utils::ElementType e)
    : InitializationException("Parameters cannot be found for the following element: " + Utils::ElementInfo::symbol(e)) {
  }
};

class ParametersDoNotExistForElementPairException : public InitializationException {
 public:
  explicit ParametersDoNotExistForElementPairException(Utils::ElementType e1, Utils::ElementType e2)
    : InitializationException("Parameters cannot be found for the following element pair: " + Utils::ElementInfo::symbol(e1) +
                              "-" + Utils::ElementInfo::symbol(e2)) {
  }
};

} // namespace Methods

} // namespace Utils
} // namespace Scine
#endif // UTILS_METHODEXCEPTIONS_H
