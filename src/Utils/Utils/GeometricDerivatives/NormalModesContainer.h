/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_NORMALMODESCONTAINER_H
#define UTILS_NORMALMODESCONTAINER_H

#include "Utils/GeometricDerivatives/NormalMode.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/MolecularTrajectory.h"
#include <vector>

namespace Scine {
namespace Utils {
class AtomCollection;
class MolecularTrajectory;

/**
 * @class NormalModesContainer NormalModesContainer.h
 * @brief Class that holds the and manages the normal modes.
 */
class NormalModesContainer {
 public:
  /**
   * @brief Returns the number of vibrational modes.
   * @return int Number of modes.
   */
  int size() const;
  /**
   * @brief Adds a normal mode to this container.
   * @param mode The mode to be added.
   */
  void add(NormalMode mode);
  /**
   * @brief Returns a const reference of the vibrational mode with index modeIndex.
   * @param modeIndex The index of the mode required.
   * @return The mode with index modeIndex.
   */
  const DisplacementCollection& getMode(int modeIndex) const;
  /**
   * @brief This function returns a molecular trajectory corresponding to a certain vibrational mode for visualization
   * purposes.
   * @param modeIndex The index of the mode required.
   * @param structure The AtomCollection of the molecular structure of interest.
   * @param scalingFactor The scaling factor applied to the mode to obtain the maximum displacement.
   * @return MolecularTrajectory The molecular trajectory representing the mode.
   */
  MolecularTrajectory getModeAsMolecularTrajectory(int modeIndex, const Utils::AtomCollection& structure,
                                                   double scalingFactor) const;
  /**
   * @brief Getter for the wave numbers corresponding to the vibrational modes.
   * @return std::vector<double> Vector of all wave numbers.
   */
  std::vector<double> getWaveNumbers() const;

 private:
  // The modes
  std::vector<NormalMode> modes_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_NORMALMODESCONTAINER_H
