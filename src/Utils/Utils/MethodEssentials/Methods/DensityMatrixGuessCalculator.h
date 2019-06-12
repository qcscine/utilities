/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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
  /**
   * @brief This function sets the number of electrons.
   * This is used in the initialize() function of the calculation method.
   * Necessary if the structure is changed, for example.
   * Previously the number of electrons was a const-reference to an unsigned int.
   * This caused some problems with gcc-7.3.0 and thus we decided to change it to
   * an int and to provide a setter method.
   * @param nElectrons new number of electrons.
   */
  virtual void setNElectrons(int nElectrons) = 0;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_DENSITYMATRIXGUESSCALCULATOR_H
