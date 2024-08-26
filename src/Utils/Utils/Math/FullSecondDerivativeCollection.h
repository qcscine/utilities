/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MATH_FULLSECONDDERIVATIVECOLLECTION_H
#define UTILS_MATH_FULLSECONDDERIVATIVECOLLECTION_H

#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Utils/Math/AutomaticDifferentiation/Second3D.h>
#include <Utils/Math/DerivativeCollection.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @class FullSecondDerivativeCollection FullSecondDerivativeCollection.h
 * @brief Container class for full second derivatives (Hessian matrix + first derivatives).
 */
class FullSecondDerivativeCollection : public DerivativeCollection {
 public:
  /**
   * @brief Constructor.
   * @param N Size of the container (number of atoms).
   */
  explicit FullSecondDerivativeCollection(int N = 0);
  /**
   * @brief Generate a GradientCollection from the full second derivatives and the displacements.
   */
  GradientCollection generateGradients(const DisplacementCollection& displacements) const;
  /**
   * @brief Generate a Gradient from the full second derivatives and the displacements.
   */
  Gradient generateGradient(int idx, const DisplacementCollection& displacements) const;
};

// Constructor with the number of atoms given as an argument
inline FullSecondDerivativeCollection::FullSecondDerivativeCollection(int N) : DerivativeCollection(N, 2) {
}

// This function generates gradients from displacements
inline GradientCollection FullSecondDerivativeCollection::generateGradients(const DisplacementCollection& displacements) const {
  GradientCollection gc(this->getNAtoms(), 3);
  const auto& refGrads = this->getReferenceGradients();
  for (int i = 0; i < this->getNAtoms(); ++i) {
    Eigen::Vector3d v = refGrads.row(i);
    for (int j = 0; j < this->getNAtoms(); ++j) {
      v += getHessianMatrix().block<3, 3>(3 * i, 3 * j) * displacements.row(j).transpose();
    }
    gc.row(i) = v;
  }

  return gc;
}

inline Gradient FullSecondDerivativeCollection::generateGradient(int /*idx*/, const DisplacementCollection& /*unused*/) const {
  throw std::runtime_error("Not implemented yet");
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_MATH_FULLSECONDDERIVATIVECOLLECTION_H
