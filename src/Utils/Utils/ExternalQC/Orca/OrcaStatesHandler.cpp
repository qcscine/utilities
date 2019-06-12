/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "OrcaStatesHandler.h"
#include "OrcaCalculator.h"
#include "OrcaState.h"
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/IO/NativeFilenames.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {

OrcaStatesHandler::OrcaStatesHandler(OrcaCalculator& calculator) : calculator_(calculator) {
}

OrcaStatesHandler::~OrcaStatesHandler() = default;

void OrcaStatesHandler::store(Utils::StateSize size) {
  auto orcaState = std::dynamic_pointer_cast<OrcaState>(getCurrentState(size));
  states_.emplace_back(orcaState);
  copyBackupFile(calculator_.getFileNameBase(), orcaState->getStringState(""));
}

void OrcaStatesHandler::load(std::shared_ptr<Utils::State> state) {
  if (!state) {
    throw EmptyStateException();
  }
  try {
    std::shared_ptr<OrcaState> orcaState = std::dynamic_pointer_cast<OrcaState>(state);
    copyBackupFile(orcaState->getStringState(""), calculator_.getFileNameBase());
  }
  catch (std::bad_cast& e) {
    throw StateNotCompatibleException();
  }
}

std::shared_ptr<Utils::State> OrcaStatesHandler::getCurrentState(Utils::StateSize size) const {
  OrcaState state(size);
  state.initialize();
  return std::make_shared<OrcaState>(state);
}

void OrcaStatesHandler::copyBackupFile(const std::string& from, const std::string& to) {
  const auto& workingDirectory = calculator_.getCalculationDirectory();
  auto fromFile = NativeFilenames::combinePathSegments(workingDirectory, from + ".gbw");
  auto toFile = NativeFilenames::combinePathSegments(workingDirectory, to + ".gbw");

  try {
    FilesystemHelpers::copyFile(fromFile, toFile);
  }
  catch (std::runtime_error& e) {
    throw StateSavingException();
  }
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
