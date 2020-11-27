/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ORCASTATE_H
#define UTILS_ORCASTATE_H

#include "Utils/IO/NativeFilenames.h"
#include "Utils/Technical/UniqueIdentifier.h"
#include <Core/BaseClasses/StateHandableObject.h>
#include <stdio.h>
#include <exception>
#include <utility>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @brief Exception for the case that a state is requested which is empty.
 */
class EmptyStateException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "ORCA state is empty and cannot be loaded.";
  }
};

/**
 * @brief Exception for a failure to save or load an ORCA state.
 */
class StateSavingException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "Failure while saving or loading ORCA state.";
  }
};

/**
 * @brief Definition of a calculation state for ORCA calculations.
 *        The calculation state is defined as a unique identifier. Only a string state is saved here.
 */
struct OrcaState final : public Core::State {
  /**
   * @brief Constructor, calls the base class constructor to initialize the size of the state.
   */
  explicit OrcaState(std::string dir) : directory(dir) {
    UniqueIdentifier id;
    stateIdentifier = id.getStringRepresentation();
  }
  /// @brief Destructor, deletes the .gbw file.
  ~OrcaState() final {
    auto file = NativeFilenames::combinePathSegments(directory, stateIdentifier + ".gbw");
    std::remove(file.c_str());
  }

  /**
   * @brief Initializer for the state.
   */
  void initialize();

  std::string directory;
  std::string stateIdentifier;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_ORCASTATE_H
