/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_IOEXCEPTIONS_H_
#define UTILS_IOEXCEPTIONS_H_

#include <stdexcept>

namespace Scine {
namespace Utils {

/**
 * @class IOException IOExceptions.h
 * @brief Base class for exceptions related to input-output.
 */
class IOException : public std::runtime_error {
 public:
  explicit IOException(const std::string& s) : std::runtime_error(s) {
  }
};
/**
 * @class InvalidFile IOExceptions.h
 * @brief  Exception when file cannot be opened.
 */
class InvalidFile : public IOException {
 public:
  explicit InvalidFile(const std::string& filename)
    : IOException("Problem when opening / creating file \"" + filename + "\"") {
  }
};
/**
 * @class InvalidAtomCountSpecification IOExceptions.h
 * @brief Exception thrown when the first line of a xyz input does not contain the number of atoms.
 */
class InvalidAtomCountSpecification : public IOException {
 public:
  explicit InvalidAtomCountSpecification(const std::string& firstLine)
    : IOException("\"" + firstLine + "\" is not a valid specification for the number of atoms") {
  }
};
/**
 * @class InvalidAtomSpecification IOExceptions.h
 * @brief Exception thrown when a line for an atom of a xyz input is invalid.
 */
class InvalidAtomSpecification : public IOException {
 public:
  explicit InvalidAtomSpecification(const std::string& atomLine)
    : IOException("\"" + atomLine + "\" is not a valid XYZ line for an atom") {
  }
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_IOEXCEPTIONS_H_
