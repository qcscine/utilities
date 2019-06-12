/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MethodInitializer.h"
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <utility>

namespace Scine {
namespace Utils {

std::string MethodInitializer::rootParameterFolder_;

void MethodInitializer::setRootParameterFolder(std::string path) {
  rootParameterFolder_ = std::move(path);
}

std::string MethodInitializer::getRootParameterFolder() {
  return rootParameterFolder_;
}

void MethodInitializer::initializeMethod(const Utils::AtomCollection& structure) {
  const auto& positions = structure.getPositions();
  const auto& elements = structure.getElements();
  initializeMethod(positions, elements);
}

void MethodInitializer::initializeMethod(const Utils::ElementTypeCollection& elementTypes) {
  initialize(elementTypes);
}

void MethodInitializer::initializeMethod(const Utils::PositionCollection& positions,
                                         const Utils::ElementTypeCollection& elementTypes) {
  initialize(positions, elementTypes);
}
} // namespace Utils
} // namespace Scine
