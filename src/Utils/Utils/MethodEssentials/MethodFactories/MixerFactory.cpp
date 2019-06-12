/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MixerFactory.h"
#include <Utils/MethodEssentials/Mixing/EdiisDiisModifier.h>
#include <Utils/MethodEssentials/Mixing/EdiisModifier.h>
#include <Utils/MethodEssentials/Mixing/FockDiisModifier.h>
#include <Utils/MethodEssentials/Mixing/charge_simple.h>
#include <Utils/MethodEssentials/Mixing/fock_simple.h>

namespace Scine {
namespace Utils {

std::vector<MixerFactory::MixerDescriptor> MixerFactory::availableMixers;
scf_mixer_t MixerFactory::defaultMixer;

const std::vector<MixerFactory::MixerDescriptor>& MixerFactory::getAvailableMixers() {
  if (availableMixers.empty())
    setAvailableMixers();
  return availableMixers;
}

void MixerFactory::setAvailableMixers() {
  availableMixers.emplace_back(scf_mixer_t::none, "No mixer");
  availableMixers.emplace_back(scf_mixer_t::fock_diis, "Fock DIIS mixer");
  availableMixers.emplace_back(scf_mixer_t::ediis, "EDIIS mixer");
  availableMixers.emplace_back(scf_mixer_t::ediis_diis, "EDIIS + DIIS mixer");
  availableMixers.emplace_back(scf_mixer_t::charge_simple, "Simple charge mixer");
  availableMixers.emplace_back(scf_mixer_t::fock_simple, "Simple Fock mixer");
  defaultMixer = scf_mixer_t::fock_diis;
}

std::unique_ptr<SCFModifier> MixerFactory::createMixer(scf_mixer_t mixerID) {
  switch (mixerID) {
    case scf_mixer_t::fock_diis:
      return std::make_unique<FockDiisModifier>();
    case scf_mixer_t::ediis:
      return std::make_unique<EdiisModifier>();
    case scf_mixer_t::ediis_diis:
      return std::make_unique<EdiisDiisModifier>();
    case scf_mixer_t::charge_simple:
      return std::make_unique<Charge_Simple>();
    case scf_mixer_t::fock_simple:
      return std::make_unique<Fock_Simple>();
    default:
      return std::unique_ptr<SCFModifier>(); // i.e. nullptr
  }
}
} // namespace Utils
} // namespace Scine
