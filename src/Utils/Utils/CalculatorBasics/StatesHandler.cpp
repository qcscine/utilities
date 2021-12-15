/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/CalculatorBasics/StatesHandler.h"

namespace Scine {
namespace Utils {

StatesHandler::StatesHandler(std::shared_ptr<Core::StateHandableObject> object) : statesHandableObject_(object) {
}

void StatesHandler::store(std::shared_ptr<Core::State> state) {
  states_.emplace_back(std::move(state));
}

void StatesHandler::store() {
  auto instance = statesHandableObject_.lock();
  if (!instance) {
    throw NoStateHandableObjectPresent();
  }
  states_.emplace_back(instance->getState());
}

void StatesHandler::load(std::shared_ptr<Core::State> state) {
  auto instance = statesHandableObject_.lock();
  if (!instance) {
    throw NoStateHandableObjectPresent();
  }
  instance->loadState(std::move(state));
}

void StatesHandler::load(int index) {
  auto instance = statesHandableObject_.lock();
  if (!instance) {
    throw NoStateHandableObjectPresent();
  }
  instance->loadState(getState(index));
}

std::shared_ptr<Core::State> StatesHandler::getState(int index) const {
  return states_.at(index);
}

std::shared_ptr<Core::State> StatesHandler::popOldestState() {
  if (!states_.empty()) {
    std::shared_ptr<Core::State> state = std::move(states_.front());
    states_.pop_front();
    return state;
  }

  throw EmptyStatesHandlerContainer();
}

std::shared_ptr<Core::State> StatesHandler::popNewestState() {
  if (!states_.empty()) {
    std::shared_ptr<Core::State> state = std::move(states_.back());
    states_.pop_back();
    return state;
  }

  throw EmptyStatesHandlerContainer();
}

void StatesHandler::clear() {
  states_.clear();
}

int StatesHandler::size() const {
  return states_.size();
}

} // namespace Utils
} // namespace Scine
