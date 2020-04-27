/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ConvergenceAcceleratorFactory.h"
#include "ChargeSimple.h"
#include "EdiisDiisModifier.h"
#include "EdiisModifier.h"
#include "FockDiisModifier.h"
#include "FockSimple.h"

namespace Scine {
namespace Utils {

std::vector<ConvergenceAcceleratorFactory::MixerDescriptor> ConvergenceAcceleratorFactory::availableMixers;
scf_mixer_t ConvergenceAcceleratorFactory::defaultMixer;

const std::vector<ConvergenceAcceleratorFactory::MixerDescriptor>& ConvergenceAcceleratorFactory::getAvailableMixers() {
  if (availableMixers.empty())
    setAvailableMixers();
  return availableMixers;
}

void ConvergenceAcceleratorFactory::setAvailableMixers() {
  availableMixers.emplace_back(scf_mixer_t::none, "No mixer");
  availableMixers.emplace_back(scf_mixer_t::fock_diis, "Fock DIIS mixer");
  availableMixers.emplace_back(scf_mixer_t::ediis, "EDIIS mixer");
  availableMixers.emplace_back(scf_mixer_t::ediis_diis, "EDIIS + DIIS mixer");
  availableMixers.emplace_back(scf_mixer_t::charge_simple, "Simple charge mixer");
  availableMixers.emplace_back(scf_mixer_t::fock_simple, "Simple Fock mixer");
  defaultMixer = scf_mixer_t::fock_diis;
}

std::unique_ptr<ScfModifier> ConvergenceAcceleratorFactory::createMixer(scf_mixer_t mixerID) {
  switch (mixerID) {
    case scf_mixer_t::fock_diis:
      return std::make_unique<FockDiisModifier>();
    case scf_mixer_t::ediis:
      return std::make_unique<EdiisModifier>();
    case scf_mixer_t::ediis_diis:
      return std::make_unique<EdiisDiisModifier>();
    case scf_mixer_t::charge_simple:
      return std::make_unique<ChargeSimple>();
    case scf_mixer_t::fock_simple:
      return std::make_unique<FockSimple>();
    default:
      return std::unique_ptr<ScfModifier>(); // i.e. nullptr
  }
}
} // namespace Utils
} // namespace Scine
