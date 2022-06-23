/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/GeometricDerivatives/NormalModeAnalysis.h"
#include "Utils/Constants.h"
#include "Utils/GeometricDerivatives/HessianUtilities.h"
#include "Utils/GeometricDerivatives/NormalMode.h"
#include "Utils/GeometricDerivatives/NormalModesContainer.h"
#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {
namespace NormalModeAnalysis {

inline NormalModesContainer calculate(HessianUtilities& diagonalizer, int nAtoms, bool normalize) {
  Eigen::VectorXd eigenvalues = diagonalizer.getInternalEigenvalues();
  Eigen::MatrixXd cartesianDisplacements = diagonalizer.getBackTransformedInternalEigenvectors(normalize);

  NormalModesContainer modesContainer;
  DisplacementCollection dc(nAtoms, 3);
  for (int i = 0; i < cartesianDisplacements.cols(); ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      dc.row(j) = cartesianDisplacements.col(i).segment(3 * j, 3);
    }

    double freq = getWaveNumber(eigenvalues[i]);
    NormalMode m(freq, dc);
    modesContainer.add(std::move(m));
  }
  return modesContainer;
}

NormalModesContainer calculateNormalModes(const HessianMatrix& hessian, const AtomCollection& atoms) {
  return calculateNormalModes(hessian, atoms.getElements(), atoms.getPositions(), true);
}

NormalModesContainer calculateNormalModes(const HessianMatrix& hessian, const ElementTypeCollection& elements,
                                          const PositionCollection& positions, bool normalize) {
  int nAtoms = elements.size();

  HessianUtilities diagonalizer(hessian, elements, positions, true);

  return calculate(diagonalizer, nAtoms, normalize);
}

NormalModesContainer calculateOrthogonalNormalModes(const HessianMatrix& hessian, const ElementTypeCollection& elements,
                                                    const PositionCollection& positions, const GradientCollection& gradient) {
  assert(gradient.size() == hessian.rows() && "Gradient dimension and hessian dimension do not match! (must be 3*N)");
  int nAtoms = elements.size();

  HessianUtilities diagonalizer(hessian, elements, positions, gradient, true);

  return calculate(diagonalizer, nAtoms, true);
}

double getWaveNumber(double value) {
  double f1 = Constants::invCentimeter_per_hartree * sqrt(Constants::u_per_electronRestMass); // for conversion to cm^-1
  double f2 = value < 0 ? -1 : 1; // for imaginary frequencies, will be shown as negative frequencies

  return f1 * f2 * std::sqrt(f2 * value);
}

} // namespace NormalModeAnalysis
} // namespace Utils
} // namespace Scine
