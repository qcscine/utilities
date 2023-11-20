/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_DENSITYMATRIXGUESSCALCULATOR_H
#define UTILS_DENSITYMATRIXGUESSCALCULATOR_H

namespace Scine {
namespace Utils {

class DensityMatrix;

/*!
 * Interface for the calculation of the density matrix guess in SCF calculations.
 */
class DensityMatrixGuessCalculator {
 public:
  virtual ~DensityMatrixGuessCalculator() = default;

  virtual DensityMatrix calculateGuess() const = 0;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_DENSITYMATRIXGUESSCALCULATOR_H
