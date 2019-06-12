/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ScfConvergenceAccelerator.h"
#include "SCFMethod.h"
#include <Utils/MethodEssentials/MethodFactories/MixerFactory.h>

namespace Scine {
namespace Utils {

ScfConvergenceAccelerator::ScfConvergenceAccelerator(SCFMethod& method)
  : method_(method), currentScheme_(scf_mixer_t::none) {
  setScfMixer(scf_mixer_t::fock_diis);
}

ScfConvergenceAccelerator::~ScfConvergenceAccelerator() {
  removeCurrentMixer();
}

void ScfConvergenceAccelerator::setScfMixer(scf_mixer_t mixer) {
  if (mixer == currentScheme_)
    return;

  removeCurrentMixer();
  setMixer(mixer);

  currentScheme_ = mixer;
}

scf_mixer_t ScfConvergenceAccelerator::getScfMixer() const {
  return currentScheme_;
}

void ScfConvergenceAccelerator::removeCurrentMixer() {
  if (activeScfMixer_)
    method_.removeModifier(activeScfMixer_);
  activeScfMixer_.reset();
}

void ScfConvergenceAccelerator::setMixer(scf_mixer_t mixer) {
  if (mixer != scf_mixer_t::none) {
    activeScfMixer_ = MixerFactory::createMixer(mixer);
    method_.addModifier(activeScfMixer_, 10);
  }
}
} // namespace Utils
} // namespace Scine
