/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_NORMALMODEANALYSIS_H
#define UTILS_NORMALMODEANALYSIS_H

#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Typenames.h>

namespace Scine {
namespace Utils {
class NormalModesContainer;

namespace NormalModeAnalysis {
/**
 * @brief Calculate the normal modes container of structure with given hessian
 *
 * @param hessian   The hessian (non-mass weighted, in cartesian coordinates).
 * @param atoms     The AtomCollection of the underlying structure.
 * @return The normal modes obtained from the mass weighted Hessian back-scaled to Cartesian coordinates summarized in a
 * container.
 */
NormalModesContainer calculateNormalModes(const HessianMatrix& hessian, const AtomCollection& atoms);

/**
 * @brief Calculate the normal modes container of structure with given hessian
 *
 * @param hessian   The hessian (non-mass weighted, in cartesian coordinates).
 * @param elements  The elements of the underlying structure.
 * @param positions The atom positions of the underlying structure.
 * @param normalize If true (default value), returns the normalized MW normal modes.
 *                  Returns the MW normal modes otherwise, with the norm given by the reduced mass.
 * @return The normal modes obtained from the mass weighted Hessian back-scaled to Cartesian coordinates summarized in a
 * container.
 */
NormalModesContainer calculateNormalModes(const HessianMatrix& hessian, const ElementTypeCollection& elements,
                                          const PositionCollection& positions, bool normalize = true);

/**
 * @brief Calculate the unnormalized normal modes of a structure with a given Hessian.
 *
 * @param hessian   The Hessian (non-mass weighted, in Cartesian coordinates).
 * @param elements  The elements of the underlying structure.
 * @param positions The atom positions of the underlying structure.
 * @return The normal modes obtained from the mass-weighted Hessian back-scaled to Cartesian coordinates summarized in a
 * container.
 */
NormalModesContainer calculateUnnormalizedNormalModes(const HessianMatrix& hessian, const ElementTypeCollection& elements,
                                                      const PositionCollection& positions);

/**
 * @brief Calculate the normal modes container of structure with given hessian s.t.
 *        it is orthogonal to a gradient.
 *
 * @param hessian   The hessian (non-mass weighted, in cartesian coordinates).
 * @param elements  The elements of the underlying structure.
 * @param positions The atom positions of the underlying structure.
 * @param gradient  The gradient to which the hessian must be orthogonal to.
 * @return The mass weighted normalmodes summerized in a container.
 */
NormalModesContainer calculateOrthogonalNormalModes(const HessianMatrix& hessian, const ElementTypeCollection& elements,
                                                    const PositionCollection& positions, const GradientCollection& gradient);
/**
 * @brief Convert internal eigenvalues of hessian to wave number in [cm^-1]
 *
 * @param value Internal eigenvalue obtained by the diagonalizer of the hessian.
 * @return Corresponding wavenumber in [cm^-1]. Imaginary wavenumbers are returned as negative numbers.
 */
double getWaveNumber(double value);

/**
 * @brief Calculate n-th harmonic inversion point for a given wave number in [cm^-1]
 *
 * @param wavenumber Wavenumber in [cm^-1]. Imaginary wavenumbers are given as negative numbers.
 * @param n Value specifying the n-th harmonic inversion point to be calculated.
 *          Usually is a positive integer, but can in principle be any positive real number.
 * @return Displacement corresponding to the n-th harmonic inversion point.
 * The displacement value should be used to displace along vibrational modes obtained from the normal mode analysis as
 * implemented above.
 */
double calculateHarmonicInversionPoint(double wavenumber, double n);

} // namespace NormalModeAnalysis
} // namespace Utils
} // namespace Scine

#endif // UTILS_NORMALMODEANALYSIS_H
