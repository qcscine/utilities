/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/CalculatorBasics/StatesHandler.h"

namespace Scine {
namespace Utils {

void StatesHandler::store(std::shared_ptr<State> state) {
  states_.emplace_back(std::move(state));
}

std::shared_ptr<State> StatesHandler::getState(int index) {
  return states_.at(index);
}

std::shared_ptr<State> StatesHandler::popOldestState() {
  if (!states_.empty()) {
    std::shared_ptr<State> state = std::move(states_.front());
    states_.pop_front();
    return state;
  }
  else {
    throw EmptyStatesHandlerContainer();
  }
}

std::shared_ptr<State> StatesHandler::popNewestState() {
  if (!states_.empty()) {
    std::shared_ptr<State> state = std::move(states_.back());
    states_.pop_back();
    return state;
  }
  else {
    throw EmptyStatesHandlerContainer();
  }
}

void StatesHandler::clear() {
  states_.clear();
}

int StatesHandler::size() const {
  return states_.size();
}

} // namespace Utils
} // namespace Scine
