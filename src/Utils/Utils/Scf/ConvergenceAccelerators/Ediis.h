/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_EDIIS_H
#define UTILS_EDIIS_H

#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {

/*!
 * Class performing the calculation of a Fock matrix based on the EDIIS algorithm.
 */
class Ediis {
 public:
  Ediis();
  void setSubspaceSize(int n);
  void setNAOs(int n);
  void setUnrestricted(bool b);
  void addMatrices(double energy, const SpinAdaptedMatrix& F, const DensityMatrix& P);
  void restart();
  SpinAdaptedMatrix getMixedFockMatrix();

 private:
  void resizeMembers();
  void updateBMatrix();

  bool unrestricted_ = false;
  int subspaceSize_;
  int nAOs_ = 0;

  int index_;
  int lastAdded_;
  int iterationNo_;

  std::vector<SpinAdaptedMatrix> fockMatrices;
  std::vector<DensityMatrix> densityMatrices;
  std::vector<double> energies;

  Eigen::MatrixXd B;
  double getBMatrixElement(int i, int j) const;
  SpinAdaptedMatrix calculateLinearCombination(const Eigen::VectorXd& coefs);
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_EDIIS_H