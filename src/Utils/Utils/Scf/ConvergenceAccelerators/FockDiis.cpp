/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "FockDiis.h"
#include <Eigen/QR>
#include <algorithm>

namespace Scine {
namespace Utils {

FockDiis::FockDiis() {
  setSubspaceSize(5);
}

void FockDiis::setUnrestricted(bool b) {
  unrestricted_ = b;
  diisError_.setUnrestricted(b);
}

void FockDiis::setSubspaceSize(int n) {
  bool resizeNeeded = n != subspaceSize_;

  subspaceSize_ = n;

  if (resizeNeeded)
    resizeMembers();
}

void FockDiis::setNAOs(int n) {
  bool resizeNeeded = n != nAOs_;

  nAOs_ = n;

  if (resizeNeeded)
    resizeMembers();
}

void FockDiis::resizeMembers() {
  fockMatrices.resize(subspaceSize_);
  diisError_.resize(subspaceSize_);
  diisStepErrors_.resize(subspaceSize_);

  overlap = Eigen::MatrixXd::Zero(nAOs_, nAOs_);

  B = Eigen::MatrixXd::Ones(subspaceSize_ + 1, subspaceSize_ + 1) * (-1);
  B(0, 0) = 0;

  rhs = Eigen::VectorXd::Zero(subspaceSize_ + 1);
  rhs(0) = -1;

  restart();
}

void FockDiis::setOverlapMatrix(const Eigen::MatrixXd& S) {
  overlap = S.selfadjointView<Eigen::Lower>();
  restart();
}

void FockDiis::restart() {
  C = Eigen::VectorXd::Zero(subspaceSize_ + 1);
  iterationNo_ = 0;
  index_ = 0;
}

void FockDiis::addMatrices(const SpinAdaptedMatrix& F, const DensityMatrix& P) {
  iterationNo_++;
  lastAdded_ = index_;

  fockMatrices[index_] = F;

  diisError_.setErrorFromMatrices(index_, F, P, overlap);
  diisStepErrors_[index_] = std::sqrt(diisError_.getError(index_, index_)) / nAOs_;

  updateBMatrix();

  index_ = (index_ + 1) % subspaceSize_;
}

void FockDiis::updateBMatrix() {
  int activeSize = iterationNo_ > subspaceSize_ ? subspaceSize_ : iterationNo_;

  // Bii element
  B(lastAdded_ + 1, lastAdded_ + 1) = diisError_.getError(lastAdded_, lastAdded_);

  // Bij elements
  for (int i = 1; i < activeSize + 1; i++) {
    if (i == lastAdded_ + 1)
      continue;
    B(lastAdded_ + 1, i) = diisError_.getError(lastAdded_, i - 1);
    B(i, lastAdded_ + 1) = B(lastAdded_ + 1, i);
  }
}

SpinAdaptedMatrix FockDiis::getMixedFockMatrix() {
  if (iterationNo_ > subspaceSize_)
    iterationNo_ = subspaceSize_;

  // If we have only one Fock matrix
  if (iterationNo_ < 2) {
    return fockMatrices[0];
  }

  C.head(iterationNo_ + 1) =
      B.block(0, 0, iterationNo_ + 1, iterationNo_ + 1).colPivHouseholderQr().solve(rhs.head(iterationNo_ + 1));

  return calculateLinearCombination();
}

SpinAdaptedMatrix FockDiis::calculateLinearCombination() {
  if (unrestricted_) {
    Eigen::MatrixXd FAlpha = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
    Eigen::MatrixXd FBeta = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
    for (int i = 0; i < iterationNo_; i++) {
      FAlpha += C(i + 1) * fockMatrices[i].alphaMatrix();
      FBeta += C(i + 1) * fockMatrices[i].betaMatrix();
    }
    return SpinAdaptedMatrix::createUnrestricted(std::move(FAlpha), std::move(FBeta));
  }
  else {
    Eigen::MatrixXd Fsol = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
    for (int i = 0; i < iterationNo_; i++)
      Fsol += C(i + 1) * fockMatrices[i].restrictedMatrix();
    return SpinAdaptedMatrix::createRestricted(std::move(Fsol));
  }
}

double FockDiis::getMaxError() const {
  int activeSize = iterationNo_ > subspaceSize_ ? subspaceSize_ : iterationNo_;
  auto maxIter = std::max_element(diisStepErrors_.begin(), diisStepErrors_.begin() + activeSize);
  return *maxIter;
}

double FockDiis::getMinError() const {
  int activeSize = iterationNo_ > subspaceSize_ ? subspaceSize_ : iterationNo_;
  auto minIter = std::min_element(diisStepErrors_.begin(), diisStepErrors_.begin() + activeSize);
  return *minIter;
}

double FockDiis::getLastError() const {
  return diisStepErrors_[lastAdded_];
}
} // namespace Utils
} // namespace Scine
