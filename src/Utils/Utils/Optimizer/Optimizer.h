/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_OPTIMIZER_H_
#define UTILS_OPTIMIZER_H_

#include "Utils/Settings.h"
#include <Eigen/Core>
#include <deque>

namespace Scine {
namespace Utils {
/**
 * @brief The base class for all Optimizers.
 *
 * Given the fact that optimizers can function in many different ways,
 * this base class is somewhat slim and only unifies the usage of auxiliary functions.\n
 * \n
 * In general the optimizers are intended to be used for internal optimizations
 * where a fixed optimizer and a fixed set of options for that optimizer are used, and
 * also those cases where user settings at runtime shall be used.\n
 * To this end, the optimizers should have their options given as public member
 * variables which can alternatively be set using the Optimizer::applySettings function.
 * In order to then allow the addition of the encoded options into a Settings object
 * a Optimizer::addSettingsDescriptors() function has to be generated in any derived
 * class.\n
 * Furthermore all optimizers shall track the number of cycles they run and also
 * call the Optimizer::triggerObservers() function once in each iteration.\n
 * \n
 * The specific function handling the optimization is dependent on the particular
 * type of optimizer, however it should have the following properties:\n
 * - it should return the final cycle count upon completion\n
 * - it should accept the (in this order):\n
 *   1. the initial parameters\n
 *   2. a template lambda function for the update of needed variables\n
 *   3. an object derived from ConvergenceCheck to signal the end of the optimization\n
 */
class Optimizer {
 public:
  /// @brief The Optimizer base class can be default constructed.
  Optimizer() = default;
  /**
   * @brief Add an observer function that will be triggered in each iteration.
   *
   * @param function A function to be executed in every loop of the optimization.
   *                 The function will have access to the current cycle count,
   *                 the current value and to a const reference of the current
   *                 parameters.
   */
  void addObserver(std::function<void(const int&, const double&, const Eigen::VectorXd&)> function) {
    _observers.push_back(function);
  }
  /**
   * @brief Clear all existing observer functions.
   *
   * For optimization problems with very fast evaluations of the underlying function
   * the removal of all observers can increase performance as the observers are given
   * as std::functions and can not be added via templates.
   */
  void clearObservers() {
    _observers.clear();
  }
  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include these optimizer options.
   *
   * This function has to be implemented in derived classes in order to allow for
   * all optimizers to be used and configured at runtime by end users.
   *
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  virtual void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const = 0;
  /**
   * @brief Updates the optimizers options with those values given in the Settings.
   *
   * This function has to be implemented in derived classes in order to allow for
   * all optimizers to be used and configured at runtime by end users.
   *
   * Note: Additional settings that are not used by the optimizer shall not generate errors,
   *       or warnings, they shall simply be ignored.
   *
   * @param settings The settings to update the option of this optimizer with.
   */
  virtual void applySettings(const Settings& settings) = 0;
  /**
   * @brief Prepares an optimizer for rerunning its optimize function
   *
   * This function is used to prepare an optimizer instance for rerunning its
   * optimize function. If a derived optimizer e.g. stores an optional Hessian
   * projection function or a custom Hessian initialisation, the specification
   * of this function shall allow for its removal.
   * This base implementation changes the cycle count the optimizer starts with
   * when the optimize function is called and clears the value memory for
   * oscillating correction.
   *
   * @param cycleNumber The cycle number the optimizer starts with.
   */
  virtual void prepareRestart(const int& cycleNumber) {
    _startCycle = cycleNumber;
    _initializedValueMemory = false;
    _valueMemory.clear();
  }
  /**
   * @brief checks if value has been oscillating over the last maxValueMemory steps
   *
   * @param value  value of the last step
   *
   * @return bool  if parameters are oscillating
   */
  bool isOscillating(const double value) {
    /* possibility to disable check by setting to zero, and storing less than 3 does not work */
    if (maxValueMemory < 3)
      return false;
    /* initialize deque with null, if first function call */
    if (!_initializedValueMemory) {
      for (int i = 0; i < maxValueMemory; ++i) {
        _valueMemory.push_back(nullptr);
      }
      _initializedValueMemory = true;
    }
    _valueMemory.pop_front();
    _valueMemory.push_back(std::make_shared<double>(value));
    /* fill deque until full, no comparisons done yet, therefore always return false */
    if (!_valueMemory[0])
      return false;
    /* return false as soon as the difference of values is not alternating signs */
    bool previousPositive = ((*_valueMemory[0] - *_valueMemory[1]) > 0.0);
    for (int i = 2; i < maxValueMemory; ++i) {
      bool thisPositive = ((*_valueMemory[i - 1] - *_valueMemory[i]) > 0.0);
      if (previousPositive == thisPositive)
        return false;
      previousPositive = thisPositive;
    }
    /* the sign of the value is alternating in saved range --> system is oscillating */
    return true;
  };
  /**
   * @brief subtracts half of the last taken stepvector, further corrections may be added in the individual optimizers
   *
   * @param steps       stepvector including scaling factor of the last step
   * @param parameters  parameters to be optimized
   *
   * @return void
   */
  void oscillationCorrection(const Eigen::VectorXd& steps, Eigen::VectorXd& parameters) {
    parameters -= (steps / 2);
  };

  // @brief the number of saved last values during optimization to check for oscillations, setting to < 3 corresponds to
  // turning the check off
  unsigned int maxValueMemory = 10;

 protected:
  /**
   * @brief Triggers all registered observers in sequence.
   * @param cycle      The current cycle counter.
   * @param value      The current value.
   * @param parameters The current parameters as optimized thus far.
   */
  void triggerObservers(const int& cycle, const double& value, const Eigen::VectorXd& parameters) {
    for (auto& function : _observers) {
      function(cycle, value, parameters);
    }
  }
  // The cycle number the optimize function starts counting with
  int _startCycle = 1;
  // A deque of the last values
  std::deque<std::shared_ptr<double>> _valueMemory;
  // if the _valueMemory has been initialized with null pointers
  bool _initializedValueMemory = false;

 private:
  // A vector of all registered observers.
  std::vector<std::function<void(const int&, const double&, const Eigen::VectorXd&)>> _observers;
};

/**
 * @brief Exception thrown when the parameters vector given to the optimizer is empty.
 */
class EmptyOptimizerParametersException : public std::runtime_error {
 public:
  explicit EmptyOptimizerParametersException() : std::runtime_error("Empty parameters passed to optimizer.") {
  }
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_OPTIMIZER_H_
