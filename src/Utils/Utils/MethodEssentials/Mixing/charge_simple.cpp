/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "charge_simple.h"
#include <Utils/MethodEssentials/Methods/SCFMethod.h>

namespace Scine {
namespace Utils {

void Charge_Simple::onIterationStart() {
  if (!initialized) {
    initialize();
    initialized = true;
  }

  addVector(m->getAtomicCharges());
  m->setAtomicCharges(extrapolate());
}

void Charge_Simple::initialize() {
  nAtoms_ = m->getNumberAtoms();
  qVectors = std::vector<std::vector<double>>(2, std::vector<double>(nAtoms_, 0.0));
  index_ = 0;
}

void Charge_Simple::addVector(const std::vector<double>& q) {
  for (int i = 0; i < nAtoms_; i++)
    qVectors[index_][i] = q[i];
  index_ = (index_ + 1) % 2;
}

const std::vector<double>& Charge_Simple::extrapolate() {
  // std::cout << "charges: " << std::endl;
  for (int a = 0; a < nAtoms_; a++) {
    // std::cout << (1.0-damping_) << "*" << qVectors[(index_+1)%2][a] << " + " << damping_ <<"*" <<qVectors[index_][a]
    // << " = ";
    qVectors[(index_ + 1) % 2][a] = (1.0 - damping_) * qVectors[(index_ + 1) % 2][a] + damping_ * qVectors[index_][a];
    // std::cout << qVectors[(index_+1)%2][a] << std::endl;
  }

  return qVectors[(index_ + 1) % 2];
}
} // namespace Utils
} // namespace Scine
