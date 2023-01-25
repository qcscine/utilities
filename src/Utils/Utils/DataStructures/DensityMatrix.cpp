/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DensityMatrix.h"
#include <cassert>

namespace Scine {
namespace Utils {

void DensityMatrix::setDensity(Matrix&& restrictedMatrix, int nElectrons) {
  assert(nElectrons >= 0);
  matrix_.setRestrictedMatrix(std::move(restrictedMatrix));
  alphaOccupation_ = nElectrons / 2.;
  betaOccupation_ = nElectrons / 2.;

  if (unrestricted_) {
    setAlphaAndBetaFromRestrictedDensity();
  }
}

void DensityMatrix::setDensity(Matrix&& alphaMatrix, Matrix&& betaMatrix, int nAlphaElectrons, int nBetaElectrons) {
  assert(nAlphaElectrons >= 0 && nBetaElectrons >= 0);
  assert(alphaMatrix.size() == betaMatrix.size() && "The alpha and beta density matrices do not have the same size.");
  alphaOccupation_ = nAlphaElectrons;
  betaOccupation_ = nBetaElectrons;
  matrix_.setRestrictedMatrix(alphaMatrix + betaMatrix);
  matrix_.setAlphaMatrix(std::move(alphaMatrix));
  matrix_.setBetaMatrix(std::move(betaMatrix));

  unrestricted_ = true;
}

void DensityMatrix::setAlphaAndBetaFromRestrictedDensity() {
  unrestricted_ = true;
  matrix_.setAlphaMatrix(matrix_.restrictedMatrix() / 2);
  matrix_.setBetaMatrix(matrix_.restrictedMatrix() / 2);
}

void DensityMatrix::resize(int nAOs) {
  matrix_.resize(nAOs);
  alphaOccupation_ = 0;
  betaOccupation_ = 0;
}

void DensityMatrix::setUnrestricted(bool b) {
  if (unrestricted_ == b) {
    return;
  }

  if (b) {
    unrestricted_ = true;
    setAlphaAndBetaFromRestrictedDensity();
  }
  else {
    unrestricted_ = false;
  }
}

auto DensityMatrix::addDensity(const DensityMatrix& rhs, double alpha) -> void {
  assert(unrestricted() == rhs.unrestricted());
  matrix_.restrictedMatrix() += alpha * rhs.matrix_.restrictedMatrix();
  if (unrestricted_) {
    matrix_.alphaMatrix() += alpha * rhs.matrix_.alphaMatrix();
    matrix_.betaMatrix() += alpha * rhs.matrix_.betaMatrix();
  }
}

auto DensityMatrix::addMatrixAlpha(const DensityMatrix::Matrix& rhs, double alpha) -> void {
  matrix_.alphaMatrix() += alpha * rhs;
}

auto DensityMatrix::addMatrixBeta(const DensityMatrix::Matrix& rhs, double alpha) -> void {
  matrix_.betaMatrix() += alpha * rhs;
}

auto DensityMatrix::addMatrixRestricted(const DensityMatrix::Matrix& rhs, double alpha) -> void {
  matrix_.restrictedMatrix() += alpha * rhs;
}

DensityMatrix& DensityMatrix::operator+=(const DensityMatrix& rhs) {
  assert(unrestricted() == rhs.unrestricted());
  matrix_.restrictedMatrix() += rhs.matrix_.restrictedMatrix();
  if (unrestricted_) {
    matrix_.alphaMatrix() += rhs.matrix_.alphaMatrix();
    matrix_.betaMatrix() += rhs.matrix_.betaMatrix();
  }
  alphaOccupation_ += rhs.alphaOccupation_;
  betaOccupation_ += rhs.betaOccupation_;
  return *this;
}

DensityMatrix DensityMatrix::operator+(const DensityMatrix& rhs) const {
  assert(unrestricted() == rhs.unrestricted());
  DensityMatrix d = *this;
  d += rhs;
  return d;
}

DensityMatrix& DensityMatrix::operator-=(const DensityMatrix& rhs) {
  assert(unrestricted() == rhs.unrestricted());
  matrix_.restrictedMatrix() -= rhs.matrix_.restrictedMatrix();
  if (unrestricted_) {
    matrix_.alphaMatrix() -= rhs.matrix_.alphaMatrix();
    matrix_.betaMatrix() -= rhs.matrix_.betaMatrix();
  }
  alphaOccupation_ -= rhs.alphaOccupation_;
  betaOccupation_ -= rhs.betaOccupation_;
  return *this;
}

DensityMatrix DensityMatrix::operator-(const DensityMatrix& rhs) const {
  assert(unrestricted() == rhs.unrestricted());
  DensityMatrix d = *this;
  d -= rhs;
  return d;
}

DensityMatrix& DensityMatrix::operator*=(double f) {
  matrix_.restrictedMatrix() *= f;
  if (unrestricted_) {
    matrix_.alphaMatrix() *= f;
    matrix_.betaMatrix() *= f;
  }
  alphaOccupation_ *= f;
  betaOccupation_ *= f;
  return *this;
}

DensityMatrix DensityMatrix::operator*(double f) {
  DensityMatrix d = *this;
  d *= f;
  return d;
}

} // namespace Utils
} // namespace Scine
