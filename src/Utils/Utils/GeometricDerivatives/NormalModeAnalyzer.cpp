/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/GeometricDerivatives/NormalModeAnalyzer.h"
#include "Utils/GeometricDerivatives/HessianUtilities.h"
#include "Utils/GeometricDerivatives/NormalMode.h"
#include "Utils/GeometricDerivatives/NormalModesContainer.h"

namespace Scine {
namespace Utils {

NormalModeAnalyzer::NormalModeAnalyzer(const HessianMatrix& hessian, const ElementTypeCollection& elements,
                                       const PositionCollection& positions)
  : hessian_(hessian), elements_(elements), positions_(positions), nAtoms_(0) {
}

NormalModesContainer NormalModeAnalyzer::calculateNormalModes() {
  nAtoms_ = elements_.size();

  HessianUtilities diagonalizer(hessian_, elements_, positions_, true);

  eigenvalues_ = diagonalizer.getInternalEigenvalues();
  cartesianDisplacements_ = diagonalizer.getBackTransformedInternalEigenvectors();

  NormalModesContainer modesContainer;
  DisplacementCollection dc(nAtoms_, 3);
  for (int i = 0; i < cartesianDisplacements_.cols(); ++i) {
    for (int j = 0; j < nAtoms_; ++j)
      dc.row(j) = cartesianDisplacements_.col(i).segment(3 * j, 3);

    double freq = getWaveNumber(eigenvalues_[i]);
    NormalMode m(freq, dc);
    modesContainer.add(std::move(m));
  }
  return modesContainer;
}

double NormalModeAnalyzer::getWaveNumber(double value) {
  double f1 = 5140.48686;         // for conversion to cm^-1
  double f2 = value < 0 ? -1 : 1; // for imaginary frequencies, will be shown as negative frequencies

  return f1 * f2 * std::sqrt(f2 * value);
}

} // namespace Utils
} // namespace Scine
