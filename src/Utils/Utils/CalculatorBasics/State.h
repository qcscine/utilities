/**
 * @file State.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_STATE_H
#define UTILS_STATE_H

#include <Eigen/Core>
#include <map>
#include <string>

namespace Scine {
namespace Utils {

/**
 * @brief enum class defining the desired size of the State to save.
 */
enum class StateSize { minimal, regular, extensive };

/**
 * @class State
 * @brief Base class for the implementation of a generic State. A state should be viewed as a checkpoint.
 * Implementation note:
 * If a class having an interface needs to save/load a state, then it has to have a polymorphic pointer to StatesHandler
 * in the interface and then populate it with a derived class, specific to the class.
 * The same is true with the polymorphic State class.
 */
class State {
 public:
  /**
   * Maps of key-value pairs containing possible states.
   */
  using MatrixState = std::map<std::string, Eigen::MatrixXd>;
  using StringState = std::map<std::string, std::string>;
  using DoubleState = std::map<std::string, double>;
  using IntState = std::map<std::string, int>;

  /**
   * @brief Constructor, constructs the state with a StateSize.
   */
  explicit State(StateSize size) : stateSize_(size){};
  /**
   * @brief Getter for the state size, i.e. the required level of information about a state.
   */
  StateSize getStateSize() const noexcept {
    return stateSize_;
  };
  /**
   * @brief Getter for an Eigen::MatrixXd state identified with a std::string.
   * This function returns a const reference to an Eigen::MatrixXd in order to avoid unnecessary copying of possibly
   * big matrices. Be careful, the state must not go out of scope!
   * @param matrixState The key for the state's key-value pair.
   * @return A const reference to the value of the key-value pair, an Eigen::MatrixXd.
   */
  virtual const Eigen::MatrixXd& getMatrixState(const std::string& matrixName) const = 0;
  /**
   * @brief Getter for a std::string state identified with a std::string.
   * This function returns a const reference to a std::string in order to avoid unnecessary copying of possibly
   * big strings. Be careful, the state must not go out of scope!
   * @param stringState The key for the state's key-value pair.
   * @return A const reference to the value of the key-value pair, a std::string.
   */
  virtual const std::string& getStringState(const std::string& stringState) const = 0;
  /**
   * @brief Getter for an integer state identified with a std::string.
   * @param intState The key for the state's key-value pair.
   * @return The value of the key-value pair, an integer.
   */
  virtual int getIntState(const std::string& intState) const = 0;
  /**
   * @brief Getter for a double state identified with a std::string.
   * @param doubleState The key for the state's key-value pair.
   * @return The value of the key-value pair, a double.
   */
  virtual double getDoubleState(const std::string& doubleState) const = 0;
  /**
   * @brief Initializer for the state, resets it to the initial state.
   */
  virtual void initialize() = 0;

 protected:
  StateSize stateSize_{StateSize::minimal};
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_STATE_H
