/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "NormalMode.h"

namespace Scine {
namespace Utils {

NormalMode::NormalMode(double waveNumber, DisplacementCollection mode)
  : waveNumber_(waveNumber), mode_(std::move(mode)) {
}

double NormalMode::getWaveNumber() const {
  return waveNumber_;
}

const DisplacementCollection& NormalMode::getMode() const {
  return mode_;
}

} // namespace Utils
} // namespace Scine
