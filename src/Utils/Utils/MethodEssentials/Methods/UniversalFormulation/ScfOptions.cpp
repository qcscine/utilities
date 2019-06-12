/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ScfOptions.h"
#include <Utils/MethodEssentials/MethodFactories/MixerFactory.h>
#include <Utils/MethodEssentials/Methods/SCFMethod.h>
#include <Utils/UniversalSettings/DescriptorCollection.h>
#include <Utils/UniversalSettings/ValueCollection.h>

namespace Scine {

using namespace Utils::UniversalSettings;

namespace Utils {

DescriptorCollection ScfOptions::getSettingDescriptors() {
  DescriptorCollection descriptors(description);

  OptionListDescriptor mixer(mixerDescription);
  mixer.addOption(noMixer);
  mixer.addOption(diis);
  mixer.addOption(ediis);
  mixer.addOption(ediisDiis);
  mixer.setDefaultOption(diis);

  DoubleDescriptor threshold(thresholdDescription);
  threshold.setMinimum(0);
  threshold.setDefaultValue(1e-5);

  IntDescriptor maxIterations(maxIterationDescription);
  maxIterations.setMinimum(1);
  maxIterations.setDefaultValue(100);

  descriptors.push_back(mixerName, std::move(mixer));
  descriptors.push_back(thresholdName, std::move(threshold));
  descriptors.push_back(maxIterationName, std::move(maxIterations));
  return descriptors;
}

void ScfOptions::applySettings(const Scine::Utils::UniversalSettings::ValueCollection& values, SCFMethod& method) {
  auto mixer = values.getString(mixerName);
  auto mixerType = convertToMixer(mixer);
  auto threshold = values.getDouble(thresholdName);
  auto maxIt = values.getInt(maxIterationName);

  method.setConvergenceCriteria(threshold);
  method.setMaxIterations(maxIt);
  method.setScfMixer(mixerType);
}

ValueCollection ScfOptions::getAppliedSettings(const SCFMethod& method) {
  auto threshold = method.getConvergenceThreshold();
  auto maxIt = method.getMaxNumberIterations();
  auto mixerType = method.getScfMixer();
  auto mixerString = convertToString(mixerType);

  ValueCollection values;
  values.addString(mixerName, mixerString);
  values.addDouble(thresholdName, threshold);
  values.addInt(maxIterationName, maxIt);
  return values;
}

scf_mixer_t ScfOptions::convertToMixer(const std::string& s) {
  if (s == noMixer) {
    return scf_mixer_t::none;
  }
  if (s == diis) {
    return scf_mixer_t::fock_diis;
  }
  if (s == ediis) {
    return scf_mixer_t::ediis;
  }
  if (s == ediisDiis) {
    return scf_mixer_t::ediis_diis;
  }

  throw std::runtime_error("Invalid value for" + std::string(mixerName));
}

std::string ScfOptions::convertToString(scf_mixer_t m) {
  switch (m) {
    case scf_mixer_t::none:
      return noMixer;
    case scf_mixer_t::fock_diis:
      return diis;
    case scf_mixer_t::ediis:
      return ediis;
    case scf_mixer_t::ediis_diis:
      return ediisDiis;
    default:
      throw std::runtime_error("Unknown conversion from scf_mixer_t to std::string");
  }
}
} // namespace Utils
} // namespace Scine
