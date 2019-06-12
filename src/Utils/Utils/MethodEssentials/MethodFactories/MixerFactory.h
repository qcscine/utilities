/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MIXERFACTORY_H
#define UTILS_MIXERFACTORY_H

#include <Utils/MethodEssentials/Methods/SCFModifier.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Scine {
namespace Utils {

enum class scf_mixer_t { none, fock_diis, ediis, ediis_diis, fock_simple, charge_simple };

/*
 * The MixerFactory class provides the register of the available Mixers
 * It is also an abstract base class for the concrete Mixer factories.
 */

class MixerFactory {
 public:
  struct MixerDescriptor {
    MixerDescriptor(scf_mixer_t ID, std::string name) : m_(ID), name_(std::move(name)) {
    }

    scf_mixer_t m_;
    std::string name_;
  };

  static scf_mixer_t defaultMixer;

  static const std::vector<MixerDescriptor>& getAvailableMixers();
  static std::unique_ptr<SCFModifier> createMixer(scf_mixer_t mixerID);

 private:
  static void setAvailableMixers();
  static std::vector<MixerDescriptor> availableMixers;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_MIXERFACTORY_H