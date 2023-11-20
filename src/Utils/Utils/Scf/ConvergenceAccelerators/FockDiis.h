/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_FOCKDIIS_H
#define UTILS_FOCKDIIS_H

#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/Scf/ConvergenceAccelerators/DiisError.h>
#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {

/*!
 * Class performing the calculation of a Fock matrix based on the direct inversion of the iterative subspace (DIIS)
 * algorithm.
 */
class FockDiis {
 public:
  FockDiis();
  void setSubspaceSize(int n);
  void setNAOs(int n);
  void setUnrestricted(bool b);
  void setOrthogonal(bool o) {
    diisError_.setOrthogonal(o);
  }
  void addMatrices(const SpinAdaptedMatrix& F, const DensityMatrix& P);
  void setOverlapMatrix(const Eigen::MatrixXd& S);
  void restart();
  SpinAdaptedMatrix getMixedFockMatrix();

  double getMaxError() const;
  double getMinError() const;
  double getLastError() const;

 private:
  void resizeMembers();
  void updateBMatrix();

  bool unrestricted_ = false;
  int subspaceSize_ = 5;
  int nAOs_ = 0;

  int index_;
  int lastAdded_;
  int iterationNo_;

  std::vector<SpinAdaptedMatrix> fockMatrices;
  DiisError diisError_;
  std::vector<double> diisStepErrors_;

  Eigen::MatrixXd overlap;

  Eigen::MatrixXd B;
  Eigen::VectorXd rhs;
  Eigen::VectorXd C;
  SpinAdaptedMatrix calculateLinearCombination();
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_FOCKDIIS_H
