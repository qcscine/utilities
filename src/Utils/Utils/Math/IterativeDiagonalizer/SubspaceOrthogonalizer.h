/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_SUBSPACEORTHOGONALIZER_H
#define UTILSOS_SUBSPACEORTHOGONALIZER_H

#include <Eigen/Core>
namespace Scine {
namespace Utils {

/**
 * @namespace SubspaceOrthogonalizer @file SubspaceOrthogonalizer.h
 * @brief Namespace containing methods to orthogonalize matrices and vectors.
 */
namespace SubspaceOrthogonalizer {
/**
 * @brief Orthogonalizes a block of the argument matrix in-place with a QR decomposition.
 * @param space A matrix whose columns need to be orthogonalized.
 * @param subspaceDimension The dimension of the block that needs to be orthogonalized.
 */
template<typename Derived>
inline void qrOrthogonalize(const Eigen::DenseBase<Derived>& space, int subspaceDimension) {
  using MatrixType = Eigen::Matrix<typename Derived::Scalar, -1, -1>;
  int spaceDimension = space.rows();
  // Calculate modified Gram-Schmidt QR decomposition
  MatrixType Q(spaceDimension, subspaceDimension);

  for (int j = 0; j < subspaceDimension; ++j) {
    Q.col(j) = space.col(j);
    for (int i = 0; i < j; ++i) {
      Q.col(j) -= Q.col(i) * (Q.col(i).adjoint() * Q.col(j));
    }
    Q.col(j).normalize();
  }
  const_cast<Eigen::DenseBase<Derived>&>(space) = Q;
}

/**
 * @brief Orthogonalizes a single vector against a subspace.
 * @param vector The vector to be orthogonalized.
 * @param subspace The subspace against which the vector is orthogonalized.
 * @param resultVector The vector where the result will be written.
 * Uses const_cast as suggested in
 * https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
 * to avoid intermediates.
 */
template<typename DerivedVecA, typename DerivedMat, typename ReturnType>
inline void orthogonalizeToSubspace(const Eigen::DenseBase<DerivedVecA>& vector, const Eigen::DenseBase<DerivedMat>& subspace,
                                    const Eigen::DenseBase<ReturnType>& returnVector) {
  const_cast<Eigen::DenseBase<ReturnType>&>(returnVector) = vector;
  for (int basis = 0; basis < subspace.cols(); ++basis) {
    const_cast<Eigen::DenseBase<ReturnType>&>(returnVector) -=
        subspace.col(basis).dot(vector.derived()) * subspace.col(basis) / subspace.col(basis).norm();
  }
}

/**
 * @brief Orthonormalizes a single vector against a subspace.
 * @param vector The vector to be orthonormalized.
 * @param subspace The subspace against which the vector is orthonormalized.
 * @param resultVector The vector where the result will be written.
 * Uses the orthogonalizeToSubspace() function.
 */
template<typename DerivedVecA, typename DerivedMat, typename ReturnType>
inline void orthonormalizeToSubspace(const Eigen::DenseBase<DerivedVecA>& vector, const Eigen::DenseBase<DerivedMat>& subspace,
                                     const Eigen::DenseBase<ReturnType>& returnVector) {
  orthogonalizeToSubspace(vector, subspace, returnVector);
  const_cast<Eigen::DenseBase<ReturnType>&>(returnVector).derived().normalize();
}

} // namespace SubspaceOrthogonalizer
} // namespace Utils
} // namespace Scine

#endif // UTILSOS_SUBSPACEORTHOGONALIZER_H
