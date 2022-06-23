/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_TURBOMOLESTATE_H
#define UTILS_TURBOMOLESTATE_H

#include "Utils/IO/FilesystemHelpers.h"
#include "Utils/IO/NativeFilenames.h"
#include "Utils/Technical/UniqueIdentifier.h"
#include <Core/BaseClasses/StateHandableObject.h>
#include <boost/filesystem.hpp>
#include <cstdio>
#include <exception>
#include <utility>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @brief Exception for a failure to save or load an TURBOMOLE state.
 */
class TurbomoleStateSavingException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "Failure while saving or loading TURBOMOLE state.";
  }
};

/**
 * @brief Definition of a calculation state for TURBOMOLE calculations.
 *        The calculation state is defined as a unique identifier. Only a string state is saved here.
 */
struct TurbomoleState final : public Core::State {
  /**
   * @brief Constructor, calls the base class constructor to initialize the size of the state.
   */
  explicit TurbomoleState(std::string dir) : directory(std::move(dir)) {
    UniqueIdentifier id;
    stateIdentifier = id.getStringRepresentation();
    FilesystemHelpers::createDirectories(stateIdentifier);
  }
  /// @brief Destructor, deletes all files.
  ~TurbomoleState() final {
    boost::filesystem::remove_all(stateIdentifier);
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

#endif // UTILS_TURBOMOLESTATE_H
