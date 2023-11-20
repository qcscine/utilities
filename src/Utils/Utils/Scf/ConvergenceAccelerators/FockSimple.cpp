/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "FockSimple.h"
#include <Utils/Scf/MethodInterfaces/ScfMethod.h>

namespace Scine {
namespace Utils {

void FockSimple::onFockCalculated() {
  if (!initialized_) {
    initialize();
    initialized_ = true;
  }
  addMatrices(m->getDensityMatrix().restrictedMatrix());
  m->setFockMatrix(SpinAdaptedMatrix::createRestricted(extrapolate()));
}

void FockSimple::initialize() {
  nAOs_ = m->getNumberAtomicOrbitals();
  fockMatrices_ = std::vector<Eigen::MatrixXd>(2, Eigen::MatrixXd::Zero(nAOs_, nAOs_));
  index_ = 0;
}

void FockSimple::addMatrices(const Eigen::MatrixXd& F) {
  fockMatrices_[index_] = F;
  index_ = (index_ + 1) % 2;
}

const Eigen::MatrixXd& FockSimple::extrapolate() {
  // Replace the Fock matrix of the current SCF cycle (which was added just before) by a linear combination of this
  // Fock matrix and the Fock matrix calculated in the previous SCF iteration.
  fockMatrices_[(index_ + 1) % 2] = (1.0 - damping_) * fockMatrices_[(index_ + 1) % 2] + damping_ * fockMatrices_[index_];
  return fockMatrices_[(index_ + 1) % 2];
}
} // namespace Utils
} // namespace Scine
