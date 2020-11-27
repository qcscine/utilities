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

NormalModesContainer calculateNormalModes(const HessianMatrix& hessian, const ElementTypeCollection& elements,
                                          const PositionCollection& positions) {
  int nAtoms = elements.size();

  HessianUtilities diagonalizer(hessian, elements, positions, true);

  Eigen::VectorXd eigenvalues = diagonalizer.getInternalEigenvalues();
  Eigen::MatrixXd cartesianDisplacements = diagonalizer.getBackTransformedInternalEigenvectors();

  NormalModesContainer modesContainer;
  DisplacementCollection dc(nAtoms, 3);
  for (int i = 0; i < cartesianDisplacements.cols(); ++i) {
    for (int j = 0; j < nAtoms; ++j)
      dc.row(j) = cartesianDisplacements.col(i).segment(3 * j, 3);

    double freq = getWaveNumber(eigenvalues[i]);
    NormalMode m(freq, dc);
    modesContainer.add(std::move(m));
  }
  return modesContainer;
}

double getWaveNumber(double value) {
  double f1 = Constants::invCentimeter_per_hartree * sqrt(Constants::u_per_electronRestMass); // for conversion to cm^-1
  double f2 = value < 0 ? -1 : 1; // for imaginary frequencies, will be shown as negative frequencies

  return f1 * f2 * std::sqrt(f2 * value);
}

} // namespace NormalModeAnalysis
} // namespace Utils
} // namespace Scine
