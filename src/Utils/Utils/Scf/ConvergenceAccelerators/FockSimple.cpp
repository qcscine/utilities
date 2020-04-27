/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "FockSimple.h"
#include <Utils/Scf/MethodInterfaces/ScfMethod.h>

namespace Scine {
namespace Utils {

void FockSimple::onFockCalculated() {
  if (!initialized) {
    initialize();
    initialized = true;
  }
  addMatrices(m->getDensityMatrix().restrictedMatrix());
  m->setFockMatrix(SpinAdaptedMatrix::createRestricted(extrapolate()));
}

void FockSimple::initialize() {
  nAOs_ = m->getNumberAtomicOrbitals();
  fockMatrices = std::vector<Eigen::MatrixXd>(2, Eigen::MatrixXd::Zero(nAOs_, nAOs_));
  index_ = 0;
}

void FockSimple::addMatrices(const Eigen::MatrixXd& F) {
  fockMatrices[index_] = F;
  index_ = (index_ + 1) % 2;
}

const Eigen::MatrixXd& FockSimple::extrapolate() {
  fockMatrices[(index_ + 1) % 2] = (1.0 - damping_) * fockMatrices[(index_ + 1) % 2] + damping_ * fockMatrices[index_];
  return fockMatrices[(index_ + 1) % 2];
}
} // namespace Utils
} // namespace Scine
