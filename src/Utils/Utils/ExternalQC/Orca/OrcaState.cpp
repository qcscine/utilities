/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "OrcaState.h"
#include <Utils/Technical/UniqueIdentifier.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {

OrcaState::OrcaState(Utils::StateSize size) : Utils::State(size) {
}

const std::string& OrcaState::getStringState(const std::string& /*stringState*/) const {
  if (!stateIdentifier_.empty())
    return stateIdentifier_;
  else
    throw EmptyStateException();
}

const Eigen::MatrixXd& OrcaState::getMatrixState(const std::string& matrixState) const {
  throw StateNotAvailableException(matrixState);
}

int OrcaState::getIntState(const std::string& intState) const {
  throw StateNotAvailableException(intState);
}

double OrcaState::getDoubleState(const std::string& doubleState) const {
  throw StateNotAvailableException(doubleState);
}

void OrcaState::initialize() {
  UniqueIdentifier id;
  stateIdentifier_ = id.getStringRepresentation();
}

bool OrcaState::hasState() const {
  return !stateIdentifier_.empty();
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
