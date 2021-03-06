/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_NORMALMODEANALYSIS_H
#define UTILS_NORMALMODEANALYSIS_H

#include <Utils/Typenames.h>

namespace Scine {
namespace Utils {
class NormalModesContainer;

namespace NormalModeAnalysis {
/**
 * @brief Calculate the normal modes container of structure with given hessian
 *
 * @param hessian   The hessian (non-mass weighted, in cartesian coordinates).
 * @param elements  The elements of the underlying structure.
 * @param positions The atom positions of the underlying structure.
 * @return The mass weighted normalmodes summerized in a container.
 */
NormalModesContainer calculateNormalModes(const HessianMatrix& hessian, const ElementTypeCollection& elements,
                                          const PositionCollection& positions);
/**
 * @brief Convert internal eigenvalues of hessian to wave number in [cm^-1]
 *
 * @param value Internal eigenvalue obtained by the diagonalizer of the hessian.
 * @return Corresponding wavenumber in [cm^-1]. Imaginary wavenumbers are returned as negative numbers.
 */
static double getWaveNumber(double value);

} // namespace NormalModeAnalysis
} // namespace Utils
} // namespace Scine

#endif // UTILS_NORMALMODEANALYSIS_H
