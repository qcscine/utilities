/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Ediis.h"
#include "EdiisCoefficientOptimizer.h"
#include <Eigen/QR>
#include <iostream>

namespace Scine {
namespace Utils {

Ediis::Ediis() {
  setSubspaceSize(6);
}

void Ediis::setUnrestricted(bool b) {
  unrestricted_ = b;
}

void Ediis::setSubspaceSize(int n) {
  bool resizeNeeded = n != subspaceSize_;

  subspaceSize_ = n;

  if (resizeNeeded)
    resizeMembers();
}

void Ediis::setNAOs(int n) {
  bool resizeNeeded = n != nAOs_;

  nAOs_ = n;

  if (resizeNeeded)
    resizeMembers();
}

void Ediis::resizeMembers() {
  fockMatrices.resize(subspaceSize_);
  densityMatrices.resize(subspaceSize_);
  energies.resize(subspaceSize_);

  B = Eigen::MatrixXd::Zero(subspaceSize_, subspaceSize_);

  restart();
}

void Ediis::restart() {
  iterationNo_ = 0;
  index_ = 0;
}

void Ediis::addMatrices(double energy, const SpinAdaptedMatrix& F, const DensityMatrix& P) {
  iterationNo_++;
  lastAdded_ = index_;

  fockMatrices[index_] = F;
  densityMatrices[index_] = P;
  energies[index_] = energy;

  updateBMatrix();

  index_ = (index_ + 1) % subspaceSize_;
}

void Ediis::updateBMatrix() {
  int activeSize = iterationNo_ > subspaceSize_ ? subspaceSize_ : iterationNo_;

  // Bii element
  B(lastAdded_, lastAdded_) = 0;

  // Bij elements
  for (int i = 0; i < activeSize; i++) {
    if (i == lastAdded_)
      continue;
    double v = getBMatrixElement(lastAdded_, i);
    B(lastAdded_, i) = v;
    B(i, lastAdded_) = B(lastAdded_, i);
  }
}

double Ediis::getBMatrixElement(int i, int j) const {
  if (unrestricted_) {
    double va = ((fockMatrices[i].alphaMatrix() - fockMatrices[j].alphaMatrix()).selfadjointView<Eigen::Lower>() *
                 (densityMatrices[i].alphaMatrix() - densityMatrices[j].alphaMatrix()))
                    .trace();
    double vb = ((fockMatrices[i].betaMatrix() - fockMatrices[j].betaMatrix()).selfadjointView<Eigen::Lower>() *
                 (densityMatrices[i].betaMatrix() - densityMatrices[j].betaMatrix()))
                    .trace();
    return (va + vb) / 2;
  }
  else {
    double v = ((fockMatrices[i].restrictedMatrix() - fockMatrices[j].restrictedMatrix()).selfadjointView<Eigen::Lower>() *
                (densityMatrices[i].restrictedMatrix() - densityMatrices[j].restrictedMatrix()))
                   .trace();
    // Divide v by two because density matrix is RHF
    return v / 2;
  }
}

SpinAdaptedMatrix Ediis::getMixedFockMatrix() {
  if (iterationNo_ > subspaceSize_)
    iterationNo_ = subspaceSize_;

  // If we have only one Fock matrix
  if (iterationNo_ < 2) {
    return fockMatrices[0];
  }
  else {
    EdiisCoefficientOptimizer opt(energies, B.block(0, 0, iterationNo_, iterationNo_));
    auto coefs = opt.getCoefficients();

    return calculateLinearCombination(coefs);
  }
}

SpinAdaptedMatrix Ediis::calculateLinearCombination(const Eigen::VectorXd& coefs) {
  if (unrestricted_) {
    Eigen::MatrixXd FAlpha = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
    Eigen::MatrixXd FBeta = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
    for (int i = 0; i < iterationNo_; i++) {
      FAlpha += coefs[i] * fockMatrices[i].alphaMatrix();
      FBeta += coefs[i] * fockMatrices[i].betaMatrix();
    }
    return SpinAdaptedMatrix::createUnrestricted(std::move(FAlpha), std::move(FBeta));
  }
  else {
    Eigen::MatrixXd Fsol = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
    for (int i = 0; i < iterationNo_; i++)
      Fsol += coefs[i] * fockMatrices[i].restrictedMatrix();
    return SpinAdaptedMatrix::createRestricted(std::move(Fsol));
  }
}
} // namespace Utils
} // namespace Scine
