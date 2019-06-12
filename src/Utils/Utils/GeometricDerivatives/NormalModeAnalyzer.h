/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_NORMALMODEANALYZER_H
#define UTILS_NORMALMODEANALYZER_H

#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {
class NormalModesContainer;

/**
 * @brief Class to calculate the normal modes of a molecule from its hessian matrix.
 */
class NormalModeAnalyzer {
 public:
  /**
   * @brief Construct a new NormalModeAnalyzer object.
   *
   * @param hessian   The hessian (non-mass weighted, in cartesian coordinates).
   * @param elements  The elements of the underlying structure.
   * @param positions The atom positions of the underlying structure.
   */
  NormalModeAnalyzer(const HessianMatrix& hessian, const ElementTypeCollection& elements, const PositionCollection& positions);
  /**
   * @brief Getter for the normal modes
   * @return NormalModesContainer The mass weighted normalmodes.
   */
  NormalModesContainer calculateNormalModes();

 private:
  /*! Returns the wave number in cm^-1 from a eigenvalue of the hessian matrix.
      If the eigenvalue is negative, will return a negative wave number instead of a complex number. */
  static double getWaveNumber(double value);

  const HessianMatrix& hessian_;
  const ElementTypeCollection& elements_;
  const PositionCollection& positions_;
  int nAtoms_;
  Eigen::MatrixXd cartesianDisplacements_;
  Eigen::VectorXd eigenvalues_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_NORMALMODEANALYZER_H
