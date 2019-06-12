/**
 * @file StatesHandler.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_STATESHANDLER_H
#define UTILS_STATESHANDLER_H

#include "Utils/CalculatorBasics/State.h"
#include <deque>
#include <exception>
#include <memory>

namespace Scine {
namespace Utils {

/**
 * @brief This exception is thrown if an empty states handler is popped.
 */
class EmptyStatesHandlerContainer : public std::exception {
  const char* what() const noexcept final {
    return "Tried to shrink the states handler container when no states were stored.";
  }
};

/**
 * @class StatesHandler
 * @brief Base class for the implementation of a generic saving and loading capability.
 * Implementation note:
 * If a class having an interface needs to save/load a state, then it has to populate the polymorphic pointers
 * to the State class with a derived class, specific to the class.
 */
class StatesHandler {
 public:
  using StatesContainer = std::deque<std::shared_ptr<State>>;

  StatesHandler() = default;
  virtual ~StatesHandler() = default;
  /**
   * @brief Store a state as the newest state.
   * @param state A pointer to the state to store internally.
   */
  void store(std::shared_ptr<State> state);
  /**
   * @brief Internally store the current state as the newest state.
   * @param size The required size of the state to save.
   */
  virtual void store(StateSize /*size*/){};
  /**
   * @brief Loads a state, i.e. it applies it.
   * @param state The state to be loaded.
   *        Can also be one of the states currently stored in the StatesHandler.
   */
  virtual void load(std::shared_ptr<State> /*state*/){};
  /**
   * @brief Gets the state with the index i.
   * @return A polymorphic smart pointer to the State at index i.
   */
  std::shared_ptr<State> getState(int index);
  /**
   * @brief Gets the current state without the need to load it.
   * @return A polymorphic shared pointer to the current State.
   */
  virtual std::shared_ptr<State> getCurrentState(StateSize /*size*/) const = 0;
  /**
   * @brief Eliminates from the underlying container and returns the oldest state recorded.
   */
  std::shared_ptr<State> popOldestState();
  /**
   * @brief Eliminates from the underlying container and returns the newest state recorded.
   */
  std::shared_ptr<State> popNewestState();
  /**
   * @brief Clear all the internally saved states.
   */
  void clear();
  /**
   * @brief Get the current size of the state storage.
   * @return The number of states currently saved.
   */
  int size() const;

 protected:
  StatesContainer states_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_STATESHANDLER_H
