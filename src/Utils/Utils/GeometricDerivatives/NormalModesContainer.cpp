/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/GeometricDerivatives/NormalModesContainer.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/MolecularTrajectory.h"
#include <cmath>

namespace Scine {
namespace Utils {

static constexpr double stepSizeForModeVisualization = 0.1;

int NormalModesContainer::size() const {
  return static_cast<int>(modes_.size());
}

void NormalModesContainer::add(NormalMode mode) {
  modes_.push_back(std::move(mode));
}

const DisplacementCollection& NormalModesContainer::getMode(int modeIndex) const {
  return modes_[modeIndex].getMode();
}

MolecularTrajectory NormalModesContainer::getModeAsMolecularTrajectory(int modeIndex, const Utils::AtomCollection& structure,
                                                                       double scalingFactor) const {
  // Construct a molecular trajectory
  Utils::MolecularTrajectory mt;
  mt.setElementTypes(structure.getElements());
  const auto& positions = structure.getPositions();
  // Get the required mode
  auto mode = getMode(modeIndex);

  // Calculate number of structures from scaling factor, this formula gives a reasonable value
  // of 40 structures for a step size of 0.1 and scaling factor of 2.
  auto numStructures = static_cast<int>(2 * scalingFactor / stepSizeForModeVisualization);

  // For loop to construct the molecular trajectory
  for (int i = 0; i < numStructures; i++) {
    mt.push_back(positions + mode * scalingFactor * sin(2 * M_PI * i / numStructures));
  }

  return mt;
}

std::vector<double> NormalModesContainer::getWaveNumbers() const {
  std::vector<double> waveNumbers;
  for (const auto& m : modes_)
    waveNumbers.push_back(m.getWaveNumber());
  return waveNumbers;
}

} // namespace Utils
} // namespace Scine
