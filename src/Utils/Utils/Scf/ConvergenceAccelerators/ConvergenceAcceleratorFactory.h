/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MIXERFACTORY_H
#define UTILS_MIXERFACTORY_H

#include <Utils/Scf/MethodInterfaces/ScfModifier.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Scine {
namespace Utils {

/**
 * @brief Enum class describing the available Scf convergence accelerators ("mixers").
 */
enum class scf_mixer_t {
  /// Element representing no accelerator at all.
  none,
  /// The Fock Direct Inversion of Iterative Subspace convergence accelerator (default).
  fock_diis,
  /// The Energy DIIS convergence accelerator.
  ediis,
  /// Combination of the previous two.
  ediis_diis,
  /// A simple fock extrapolation scheme. Does not work that well.
  fock_simple,
  /// Very simple charge extrapolation scheme. Does not work that well.
  charge_simple
};

/**
 * @brief The ConvergenceAcceleratorFactory class provides the register of the available convergence accelerators.
 * This class also creates pointer to the desired Scf convergence accelerators.
 * In the code convergence accelerators are often referred to as "mixers".
 */
class ConvergenceAcceleratorFactory {
 public:
  /**
   * @brief Small struct used as a map between string name and enum ID.
   */
  struct MixerDescriptor {
    MixerDescriptor(scf_mixer_t ID, std::string name) : m_(ID), name_(std::move(name)) {
    }

    scf_mixer_t m_;
    std::string name_;
  };

  /**
   * @brief Returns the enum class member corresponding to the default convergence accelerator
   */
  static scf_mixer_t defaultMixer;

  static const std::vector<MixerDescriptor>& getAvailableMixers();

  /**
   * @brief Factory method to create a convergence accelerator.
   * @param mixerID The corresponding enum class member.
   * @return A std::unique_ptr<ScfModifier> representing a polymorphic pointer to a convergence accelerator.
   */
  static std::unique_ptr<ScfModifier> createMixer(scf_mixer_t mixerID);

 private:
  static void setAvailableMixers();
  static std::vector<MixerDescriptor> availableMixers;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_MIXERFACTORY_H