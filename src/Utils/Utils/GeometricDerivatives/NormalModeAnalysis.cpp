/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/GeometricDerivatives/NormalModeAnalysis.h"
#include "Utils/Constants.h"
#include "Utils/GeometricDerivatives/HessianUtilities.h"
#include "Utils/GeometricDerivatives/NormalMode.h"
#include "Utils/GeometricDerivatives/NormalModesContainer.h"
#include "Utils/Geometry/Utilities/Transformations.h"
#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {
namespace NormalModeAnalysis {

inline NormalModesContainer calculate(HessianUtilities& diagonalizer, int nAtoms, bool normalize) {
  const Eigen::VectorXd eigenvalues = diagonalizer.getInternalEigenvalues();
  const Eigen::MatrixXd cartesianDisplacements = diagonalizer.getBackTransformedInternalEigenvectors(normalize);

  NormalModesContainer modesContainer;
  DisplacementCollection dc(nAtoms, 3);
  for (int i = 0; i < cartesianDisplacements.cols(); ++i) {
    for (long j = 0; j < nAtoms; ++j) {
      dc.row(j) = cartesianDisplacements.col(i).segment(3 * j, 3);
    }

    double freq = getWaveNumber(eigenvalues[i]);
    NormalMode m(freq, dc);
    modesContainer.add(std::move(m));
  }
  return modesContainer;
}

inline NormalModesContainer calculateFromPartial(HessianUtilities& diagonalizer, const std::vector<int>& partialIndices,
                                                 int nSuperAtoms, int nPartialAtoms, int expectedModes, bool normalize) {
  if (nPartialAtoms > nSuperAtoms) {
    throw std::runtime_error("Number of partial atoms must be smaller or equal to the number of super atoms.");
  }
  const Eigen::VectorXd eigenvalues = diagonalizer.getInternalEigenvalues();
  const Eigen::MatrixXd cartesianDisplacements = diagonalizer.getBackTransformedInternalEigenvectors(normalize);

  NormalModesContainer modesContainer;
  for (int i = 0; i < cartesianDisplacements.cols(); ++i) {
    DisplacementCollection dc = DisplacementCollection::Zero(nSuperAtoms, 3);
    for (long j = 0; j < nPartialAtoms; ++j) {
      // fill displacement values only for atoms part of the partial Hessian, otherwise keep zeros
      dc.row(partialIndices[j]) = cartesianDisplacements.col(i).segment(3 * j, 3);
    }

    const double freq = getWaveNumber(eigenvalues[i]);
    NormalMode m(freq, dc);
    modesContainer.add(std::move(m));
  }
  // add zero modes for atoms not part of the partial Hessian
  for (long i = cartesianDisplacements.cols(); i < expectedModes; ++i) {
    NormalMode m(0.0, DisplacementCollection::Zero(nSuperAtoms, 3));
    modesContainer.add(std::move(m));
  }
  return modesContainer;
}

NormalModesContainer calculateNormalModes(const HessianMatrix& hessian, const AtomCollection& atoms) {
  return calculateNormalModes(hessian, atoms.getElements(), atoms.getPositions(), true);
}

NormalModesContainer calculateNormalModes(const PartialHessian& hessian, const AtomCollection& atoms) {
  return calculateNormalModes(hessian, atoms.getElements(), atoms.getPositions(), true);
}

NormalModesContainer calculateNormalModes(const HessianMatrix& hessian, const ElementTypeCollection& elements,
                                          const PositionCollection& positions, bool normalize) {
  int nAtoms = elements.size();

  HessianUtilities diagonalizer(hessian, elements, positions, true);

  return calculate(diagonalizer, nAtoms, normalize);
}

NormalModesContainer calculateNormalModes(const PartialHessian& hessian, const ElementTypeCollection& elements,
                                          const PositionCollection& positions, bool normalize) {
  const int nSuperAtoms = elements.size();
  const int nPartialAtoms = hessian.getNumberOfAtoms();
  const auto& partialAtoms = hessian.getPartialAtoms(elements, positions);

  // determine number of modes of super system
  auto rotoTranslation = Geometry::Transformations::calculateTranslationAndRotationModes(positions, elements);
  int nSuperModes = static_cast<int>(rotoTranslation.rows());

  HessianUtilities diagonalizer(hessian.getMatrix(), partialAtoms.getElements(), partialAtoms.getPositions(), true);
  return calculateFromPartial(diagonalizer, hessian.getIndices(), nSuperAtoms, nPartialAtoms, nSuperModes, normalize);
}

NormalModesContainer calculateOrthogonalNormalModes(const HessianMatrix& hessian, const ElementTypeCollection& elements,
                                                    const PositionCollection& positions, const GradientCollection& gradient) {
  assert(gradient.size() == hessian.rows() && "Gradient dimension and Hessian dimension do not match! (must be 3*N)");
  int nAtoms = elements.size();

  HessianUtilities diagonalizer(hessian, elements, positions, gradient, true);

  return calculate(diagonalizer, nAtoms, true);
}

NormalModesContainer calculateOrthogonalNormalModes(const PartialHessian& hessian, const ElementTypeCollection& elements,
                                                    const PositionCollection& positions, const GradientCollection& gradient) {
  assert(gradient.size() == hessian.getMatrix().rows() &&
         "Gradient dimension and Hessian dimension do not match! (must be 3*N)");
  const int nSuperAtoms = elements.size();
  const int nPartialAtoms = hessian.getNumberOfAtoms();
  const auto& partialAtoms = hessian.getPartialAtoms(elements, positions);

  // determine number of modes of super system
  auto rotoTranslation = Geometry::Transformations::calculateTranslationAndRotationModes(positions, elements);
  int nSuperModes = static_cast<int>(rotoTranslation.rows());

  HessianUtilities diagonalizer(hessian.getMatrix(), partialAtoms.getElements(), partialAtoms.getPositions(), gradient, true);

  return calculateFromPartial(diagonalizer, hessian.getIndices(), nSuperAtoms, nPartialAtoms, nSuperModes, true);
}

double getWaveNumber(double value) {
  double f1 = Constants::invCentimeter_per_hartree * sqrt(Constants::u_per_electronRestMass); // for conversion to cm^-1
  double f2 = value < 0 ? -1 : 1; // for imaginary frequencies, will be shown as negative frequencies

  return f1 * f2 * std::sqrt(f2 * value);
}

double calculateHarmonicInversionPoint(double wavenumber, double n) {
  // Formula: $E_n = E_{pot}$, so $\hbar \omega (n + 1/2) = 1/2 m \omega^2 q^2$
  // Then solve for q, use mass-weighted coordinates.
  // For this calculation, everything is converted to Hartree atomic units.
  // Imaginary wavenumbers are treated as regular ones.
  double wavenumberInHartree = abs(wavenumber) * Constants::hartree_per_invCentimeter;
  double dispOfInversionPoint = sqrt((2. * (n + 0.5)) / wavenumberInHartree);
  // The mass-weighting in the normal mode analysis also needs to be converted to Hartree atomic units.
  dispOfInversionPoint /= sqrt(Constants::electronRestMass_per_u);
  return dispOfInversionPoint;
}

} // namespace NormalModeAnalysis
} // namespace Utils
} // namespace Scine
