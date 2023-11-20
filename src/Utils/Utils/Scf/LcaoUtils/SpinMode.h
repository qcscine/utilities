/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SPINMODE_H_
#define UTILS_SPINMODE_H_

#include <stdexcept>
#include <string>

namespace Scine {
namespace Utils {

/**
 * @brief Enum Class specifying the a spin mode such as 'restricted'
 *
 * Restricted represents a restricted spin
 * Unrestricted represents an unrestricted spin
 * RestrictedOpenShell represents a spin restricted open shell
 * Any represents restricted spin if multiplicity == 1 else unrestricted
 * None represents that spin is not considered at all
 */
enum class SpinMode { Restricted, Unrestricted, RestrictedOpenShell, Any, None };

class SpinModeInterpreter {
 public:
  // @brief static functions only
  SpinModeInterpreter() = delete;
  static SpinMode getSpinModeFromString(const std::string& spinMode) {
    if (spinMode == "restricted") {
      return SpinMode::Restricted;
    }
    if (spinMode == "unrestricted") {
      return SpinMode::Unrestricted;
    }
    if (spinMode == "restricted_open_shell") {
      return SpinMode::RestrictedOpenShell;
    }
    if (spinMode == "any") {
      return SpinMode::Any;
    }
    if (spinMode == "none") {
      return SpinMode::None;
    }
    throw std::logic_error("Unknown spin mode " + spinMode);
  };

  static std::string getStringFromSpinMode(const SpinMode& spinMode) {
    if (spinMode == SpinMode::Restricted) {
      return "restricted";
    }
    if (spinMode == SpinMode::Unrestricted) {
      return "unrestricted";
    }
    if (spinMode == SpinMode::RestrictedOpenShell) {
      return "restricted_open_shell";
    }
    if (spinMode == SpinMode::Any) {
      return "any";
    }
    if (spinMode == SpinMode::None) {
      return "none";
    }
    throw std::logic_error("Unknown string representation for this spin mode.");
  };
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_SPINMODE_H_
