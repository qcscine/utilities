/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SPECTROSCOPY_HESSIANUTILITIES_H
#define UTILS_SPECTROSCOPY_HESSIANUTILITIES_H

#include <Utils/Geometry/ElementTypes.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <memory>

namespace Scine {
namespace Utils {

// As an alternative to deciding for the entire class on whether to do the analyses in a mass weighted framework
// one could also compute both transformation matrices in the constructor and offer the choice when calling specific
// functions However, this is likely to result into unnecessary effort bc. always two orthogonalization procedures would
// be carried out
/**
 * @brief A utility class for Hessians allowing easier access to eigenvalues and eigenvectors
 *        of transformed versions.
 * @param hessian Hessian in cartesian coordinates (not mass-weighted, also not if massWeighted is set to true)
 * @param massWeighted If supposed to yield mass-weighted results set to true
 */
class HessianUtilities {
 public:
  HessianUtilities(const Eigen::MatrixXd& hessian, const ElementTypeCollection& elements,
                   const PositionCollection& positions, bool massWeighted);
  HessianUtilities(const Eigen::MatrixXd& hessian, const ElementTypeCollection& elements,
                   const PositionCollection& positions, const GradientCollection& gradient, bool massWeighted);
  /**
   * @brief Brief signals that the referenced Hessian has changed and deletes cached data.
   */
  void hessianUpdate();
  /**
   * @brief Replaces the current Hessian and deletes cached data.
   * @param hessian The new Hessian reference.
   */
  void hessianUpdate(const HessianMatrix& hessian);
  /**
   * @brief Get the transformation matrix removing rotational and translational contributions.
   * @return const Eigen::MatrixXd& The transformation matrix.
   */
  const Eigen::MatrixXd& getTransformationMatrix() const;
  /**
   * @brief Getter for the eigenvalues of the transformed Matrix.
   *
   * The eigenvalues are lazily evaluated and cached internally.
   *
   * @return const Eigen::VectorXd&
   */
  const Eigen::VectorXd& getInternalEigenvalues();
  /**
   * @brief Getter for the eigenvectors of the transformed Matrix.
   *
   * The eigenvectors are lazily evaluated and cached internally.
   *
   * @return const Eigen::MatrixXd&
   */
  const Eigen::MatrixXd& getInternalEigenvectors();
  /**
   * @brief Get the back-transformed internal eigenvectors without rotation and translation modes in cartesian
   * coordinates.
   *
   * @return Eigen::MatrixXd
   */
  Eigen::MatrixXd getBackTransformedInternalEigenvectors();
  /**
   * @brief Returns the transformed (pseudo-internal coordinates) Hessian.
   * @return const Eigen::MatrixXd& The Hessian.
   */
  Eigen::MatrixXd getInternalHessian() const;

 private:
  bool _massWeighted;
  std::reference_wrapper<const Eigen::MatrixXd> _hessian;
  const ElementTypeCollection& _elements;
  // Cached transformation
  Eigen::MatrixXd _transformation;
  // Calculates the lazy part of the cached data
  void calculateInternal();
  // Optional gradient member
  std::unique_ptr<GradientCollection> _gradient;
  // Internal lazy cache
  std::unique_ptr<Eigen::VectorXd> _internalEValues;
  std::unique_ptr<Eigen::MatrixXd> _internalEVectors;
};

} // namespace Utils
} // namespace Scine

class EmptyInternalHessianException : public std::runtime_error {
 public:
  explicit EmptyInternalHessianException()
    : std::runtime_error("Transformation to internals results into an empty Hessian.") {
  }
};
#endif // UTILS_SPECTROSCOPY_HESSIANUTILITIES_H
