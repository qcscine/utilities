/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_CP2KSTATE_H
#define UTILS_CP2KSTATE_H

#include "Utils/IO/NativeFilenames.h"
#include "Utils/Technical/UniqueIdentifier.h"
#include <Core/BaseClasses/StateHandableObject.h>
#include <cstdio>
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
    return "CP2K state is empty and cannot be loaded.";
  }
};

/**
 * @brief Exception for a failure to save or load an CP2K state.
 */
class Cp2kStateSavingException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "Failure while saving or loading CP2K state.";
  }
};

/**
 * @brief Definition of a calculation state for CP2K calculations.
 *        The calculation state is defined as a unique identifier. Only a string state is saved here.
 */
struct Cp2kState final : public Core::State {
  /**
   * @brief Constructor, calls the base class constructor to initialize the size of the state.
   */
  explicit Cp2kState(std::string dir) : directory(std::move(dir)) {
    UniqueIdentifier id;
    stateIdentifier = id.getStringRepresentation();
  }
  /// @brief Destructor, deletes the .gbw file.
  ~Cp2kState() final {
    auto file = NativeFilenames::combinePathSegments(directory, stateIdentifier + "-RESTART.wfn");
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

#endif // UTILS_CP2KSTATE_H
