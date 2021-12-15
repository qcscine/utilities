/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_NORMALMODE_H
#define UTILS_NORMALMODE_H

#include <Utils/Typenames.h>

namespace Scine {
namespace Utils {

/**
 * @brief Container for a single normal mode.
 */
class NormalMode {
 public:
  /**
   * @brief Construct a new Normal Mode object
   *
   * @param waveNumber The wave number, in cm^-1
   * @param mode       The mode, obtained from a mass-weighted Hessian, back-scaled to cartesian coordinates.
   */
  explicit NormalMode(double waveNumber, DisplacementCollection mode);
  /**
   * @brief Get the wave number.
   *
   * @return double  The wave number, in cm^-1.
   */
  double getWaveNumber() const;
  /**
   * @brief Getter for the Mode.
   * @return const DisplacementCollection& The mode.
   */
  const DisplacementCollection& getMode() const;

 private:
  double waveNumber_;
  DisplacementCollection mode_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_NORMALMODE_H
