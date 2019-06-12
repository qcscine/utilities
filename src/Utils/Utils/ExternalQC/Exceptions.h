/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_EXCEPTIONS_H
#define UTILS_EXTERNALQC_EXCEPTIONS_H

#include <stdexcept>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @class Exception Exceptions.h
 * @brief Base exception.
 */
class Exception : public std::runtime_error {
 public:
  /**
   * @brief Constructor.
   */
  explicit Exception(const std::string& s) : runtime_error(s) {
  }
};

/**
 * @class UnsuccessfulSystemCommand Exceptions.h
 * @brief Exception thrown when there is an error while executing the command for the external program.
 */
class UnsuccessfulSystemCommand : public Exception {
 public:
  /**
   * @brief Constructor.
   */
  explicit UnsuccessfulSystemCommand(const std::string& command, const std::string& inputFile, const std::string& outputFile)
    : Exception(createErrorMessage(command, inputFile, outputFile)) {
  }

 private:
  static std::string createErrorMessage(const std::string& command, const std::string& inputFile, const std::string& outputFile) {
    std::string message = "Error executing the command \"" + command;
    if (!inputFile.empty()) {
      message += " < " + inputFile;
    }
    if (!outputFile.empty()) {
      message += " > " + outputFile;
    }
    return message;
  }
};

/**
 * @class OutputFileParsingError Exceptions.h
 * @brief Exception thrown for errors during output file parsing.
 */
class OutputFileParsingError : public Exception {
 public:
  /**
   * @brief Constructor.
   */
  explicit OutputFileParsingError(const std::string& s) : Exception(s) {
  }
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_EXCEPTIONS_H