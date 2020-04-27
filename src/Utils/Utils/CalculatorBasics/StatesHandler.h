/**
 * @file StatesHandler.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_STATESHANDLER_H
#define UTILS_STATESHANDLER_H

#include <Core/BaseClasses/StateHandableObject.h>
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
 * @brief This exception is thrown if an empty states handler is popped.
 */
class NoStateHandableObjectPresent : public std::exception {
  const char* what() const noexcept final {
    return "The pointer in the StatesHandler is empty. Please call the constructor with a pointer to "
           "an instance of the StateHandableObject.";
  }
};
/**
 * @class StatesHandler
 * @brief Class responsible for saving, storing and loading object-specific states.
 * Core::State is an empty class. This base class is then implemented to have some
 * meaningful definition of a state for an application. I could be, for instance,
 * a density matrix for a Core::Calculator, a partially converged Markov Chain
 * Montecarlo, a checkpoint,... .
 * A Core::StateHandableObject interface exposes the functions for getting and
 * loading a single state.
 * The StatesHandler job is that of storing a list of these states and
 * freely switching between them.
 * This class is stand-alone, but can also be specialized to suit certain needs,
 * this is why it exposes a virtual destructor.
 */
class StatesHandler {
 public:
  /**
   * @brief The data structure in which states are stored.
   */
  using StatesContainer = std::deque<std::shared_ptr<Core::State>>;
  /**
   * @brief Constructor of the class taking a StateHandableObject instance.
   */
  explicit StatesHandler(std::shared_ptr<Core::StateHandableObject> object = nullptr);
  /**
   * @brief Virtual default destructor.
   */
  virtual ~StatesHandler() = default;
  /**
   * @brief Store an externally generated state as the newest state.
   * @param state A pointer to the state to store.
   */
  void store(std::shared_ptr<Core::State> state);
  /**
   * @brief Stores the current state of the Core::StateHandableObject instance in this class.
   */
  void store();
  /**
   * @brief Loads an externally generated state.
   * @param state The state to be loaded.
   * The specifics of how exactly the state is stored is not implemented in this class, but
   * is in the Core::StateHandableObject derived class.
   * This just calls the load() method in Core::StateHandableObject instance of this class.
   */
  void load(std::shared_ptr<Core::State> state);
  /**
   * @brief Loads an internally stored state.
   * @param index The index of the state to be loaded.
   * The specifics of how exactly the state is stored is not implemented in this class, but
   * is in the Core::StateHandableObject derived class.
   * This just extracts the state and calls the load() method in Core::StateHandableObject instance of this class.
   */
  void load(int index);
  /**
   * @brief Gets the state at some index in the StatesContainer.
   * @return A polymorphic smart pointer to the Core::State at some index.
   */
  std::shared_ptr<Core::State> getState(int index) const;
  /**
   * @brief Eliminates the oldest state from the StatesContainer and returns it.
   */
  std::shared_ptr<Core::State> popOldestState();
  /**
   * @brief Eliminates the newest state from the underlying container and returns it.
   */
  std::shared_ptr<Core::State> popNewestState();
  /**
   * @brief Clears all the internally saved states.
   */
  void clear();
  /**
   * @brief Gets the current size of the state storage.
   */
  int size() const;

 protected:
  std::weak_ptr<Core::StateHandableObject> statesHandableObject_;
  StatesContainer states_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_STATESHANDLER_H
