/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ORCASTATE_H
#define UTILS_ORCASTATE_H

#include <Utils/CalculatorBasics/State.h>
#include <exception>
#include <utility>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @brief Exception for the case that a state is requested which is not available.
 */
class StateNotAvailableException : public std::exception {
 public:
  explicit StateNotAvailableException(std::string stateName)
    : stateName_("ORCA state " + std::move(stateName) + " is not available.") {
  }
  const char* what() const noexcept final {
    return stateName_.c_str();
  }

 private:
  std::string stateName_;
};
class EmptyStateException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "ORCA state is empty and cannot be loaded.";
  }
};

/**
 * @brief Definition of a calculation state for ORCA calculations.
 *        The calculation state is defined as a unique identifier. Only a string state is saved here.
 */
struct OrcaState : public Utils::State {
  /**
   * @brief Constructor, calls the base class constructor to initialize the size of the state.
   */
  explicit OrcaState(Utils::StateSize size);

  /**
   * @brief Getter for a Eigen::MatrixXd state identified with a std::string.
   * @param matrixState The key for the state's key-value pair.
   * @return The value of the key-value pair, an Eigen::MatrixXd.
   */
  const Eigen::MatrixXd& getMatrixState(const std::string& matrixState) const final;
  /**
   * @brief Getter for a std::string state identified with a std::string.
   * @param stringState The key for the state's key-value pair.
   * @return The value of the key-value pair, a std::string.
   */
  const std::string& getStringState(const std::string& stringState) const final;
  /**
   * @brief Getter for a integer state identified with a std::string.
   * @param intState The key for the state's key-value pair.
   * @return The value of the key-value pair, an integer.
   */
  int getIntState(const std::string& intState) const final;
  /**
   * @brief Getter for a double state identified with a std::string.
   * @param doubleState The key for the state's key-value pair.
   * @return The value of the key-value pair, a double.
   */
  double getDoubleState(const std::string& doubleState) const final;

  /**
   * @brief Initializer for the state.
   */
  void initialize() final;

  /**
   * @brief Returns whether a state has been initialized
   */
  bool hasState() const;

 private:
  std::string stateIdentifier_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_ORCASTATE_H
