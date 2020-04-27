/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_DIISERROR_H
#define UTILS_DIISERROR_H

#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {

class SpinAdaptedMatrix;
class DensityMatrix;

/*!
 * Class to handle the error vectors used in DIIS.
 */
class DiisError {
 public:
  void setErrorFromMatrices(int index, const SpinAdaptedMatrix& fock, const DensityMatrix& density,
                            const Eigen::MatrixXd& overlap);
  void setOrthogonal(bool o);
  void setUnrestricted(bool u);
  void resize(unsigned subspaceSize);
  double getError(int i1, int i2) const;

 private:
  Eigen::MatrixXd calculateUnrestrictedErrorMatrix(const SpinAdaptedMatrix& fock, const DensityMatrix& density,
                                                   const Eigen::MatrixXd& overlap) const;
  Eigen::MatrixXd calculateRestrictedErrorMatrix(const SpinAdaptedMatrix& fock, const DensityMatrix& density,
                                                 const Eigen::MatrixXd& overlap) const;
  Eigen::MatrixXd calculateErrorMatrix(const Eigen::MatrixXd& fock, const Eigen::MatrixXd& overlap,
                                       const Eigen::MatrixXd& density) const;
  Eigen::MatrixXd calculateOrthogonalErrorMatrix(const Eigen::MatrixXd& fock, const Eigen::MatrixXd& density) const;

  bool orthogonal_ = false;
  bool unrestricted_ = false;
  std::vector<Eigen::MatrixXd> errorMatrices;
};

inline void DiisError::setOrthogonal(bool o) {
  orthogonal_ = o;
}

inline void DiisError::setUnrestricted(bool u) {
  unrestricted_ = u;
}

} // namespace Utils
} // namespace Scine
#endif // UTILS_DIISERROR_H