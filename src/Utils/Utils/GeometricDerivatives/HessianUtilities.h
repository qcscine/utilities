/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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

/**
 * @brief A utility class for Hessians allowing easier access to eigenvalues and eigenvectors
 *        of transformed versions.
 */
class HessianUtilities {
 public:
  HessianUtilities(const Eigen::MatrixXd& hessian, const ElementTypeCollection& elements, const PositionCollection& positions);
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
   * @param massWeighted If true returns the massweighted transformed eigenvalues.
   * @return const Eigen::VectorXd&
   */
  const Eigen::VectorXd& getInternalEigenvalues(bool massWeighted = false);
  /**
   * @brief Getter for the eigenvectors of the transformed Matrix.
   *
   * The eigenvectors are lazily evaluated and cached internally.
   *
   * @return const Eigen::MatrixXd&
   */
  const Eigen::MatrixXd& getInternalEigenvectors();
  /**
   * @brief Get the back-transformed internal eigenvectors without rotation and translation modes.
   *
   * @param massWeighted If true returns the massweighted backtransformed eigenvectors
   *                     (without rotational and translational degrees of freedom).
   * @return Eigen::MatrixXd
   */
  Eigen::MatrixXd getBackTransformedInternalEigenvectors(bool massWeighted);
  /**
   * @brief Returns the transformed (pseudo-internal coordinates) Hessian.
   * @param massWeighted If true returns the massweighted transformed Hessian.
   * @return const Eigen::MatrixXd& The Hessian.
   */
  Eigen::MatrixXd getInternalHessian(bool massWeighted = false) const;

 private:
  std::reference_wrapper<const Eigen::MatrixXd> _hessian;
  const ElementTypeCollection& _elements;
  const PositionCollection& _positions;
  // Cached transformation
  Eigen::MatrixXd _transformation;
  // Calculates the lazy part of the cached data
  void calculateInternal(bool massWeighted = false);
  // Internal lazy cache
  std::unique_ptr<Eigen::VectorXd> _mwinternalEValues;
  std::unique_ptr<Eigen::MatrixXd> _mwinternalEVectors;
  std::unique_ptr<Eigen::VectorXd> _internalEValues;
  std::unique_ptr<Eigen::MatrixXd> _internalEVectors;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_SPECTROSCOPY_HESSIANUTILITIES_H
