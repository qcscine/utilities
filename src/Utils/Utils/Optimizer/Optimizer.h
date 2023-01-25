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
#include <Eigen/Dense>
#include <deque>

namespace Scine {
namespace Utils {

/**
 * @brief Exception thrown when the Hessian has no negative eigenvalues.
 */
struct NoNegativeEigenValueException : public std::runtime_error {
  explicit NoNegativeEigenValueException() : std::runtime_error("No more negative eigenvalues, optimization failed.") {
  }
};

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
    _observers.push_back(std::move(function));
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
  virtual void prepareRestart(const int cycleNumber) {
    _startCycle = cycleNumber;
    _valueMemory.clear();
    _prevFollowedEigenvector.resize(0);
  }
  /**
   * @brief Reset an already initialized optimizer
   *
   * This function is used to reset an optimizer instance for rerunning its
   * optimize function. If a derived optimizer e.g. stores an optional Hessian
   * projection function or a custom Hessian initialisation, the specification
   * of this function shall allow for its removal.
   * This base implementation changes the cycle count the optimizer starts with
   * back to 1 when the optimize function is called and clears the value memory for
   * oscillating correction.
   */
  virtual void reset() {
    this->prepareRestart(1);
  };

  /**
   * @brief Read the current value of the cycle counter variable.
   *
   * @return The current cycle number.
   */
  int getCycle() const {
    return _cycle;
  };

  /**
   * @brief checks if value has been oscillating over the last maxValueMemory steps
   *
   * @param value  value of the last step
   *
   * @return bool  if parameters are oscillating
   */
  bool isOscillating(const double value) {
    /* possibility to disable check by setting to zero, and storing less than 3 does not work */
    if (maxValueMemory < 3) {
      return false;
    }
    // Pop the oldest value if the deque is at memory capacity
    if (_valueMemory.size() == maxValueMemory) {
      _valueMemory.pop_front();
    }
    _valueMemory.push_back(value);
    /* wait to judge until memory is full, therefore always return false */
    if (_valueMemory.size() < maxValueMemory) {
      return false;
    }
    /* return false as soon as the difference of values is not alternating signs */
    bool previousPositive = ((_valueMemory[0] - _valueMemory[1]) > 0.0);
    if (std::abs(_valueMemory[0] - _valueMemory[1]) < 1e-12)
      return false;
    for (unsigned i = 2; i < maxValueMemory; ++i) {
      const bool thisPositive = ((_valueMemory[i - 1] - _valueMemory[i]) > 0.0);
      if (previousPositive == thisPositive) {
        return false;
      }
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
    constrainedAdd(parameters, -steps / 2);
  };

  /**
   * @brief Determine stepvector to maximize energy along an eigenvector and minimize along all others
   *
   * @param gradients          The gradients.
   * @param hessian            The Hessian.
   * @param modeToMaximize     The number of the eigenvector that shall be maximized starting from zero with the
   *                           lowest eigenvalue. This number can be adjusted based on a previously followed eigenvector
   *                           based on cosine similarity between eigenvectors.
   * @return Eigen::VectorXd   Returns the stepvector.
   */
  Eigen::VectorXd modeMaximizationWithHessian(const Eigen::VectorXd& gradients, const Eigen::MatrixXd& hessian,
                                              int& modeToMaximize) {
    /* Decompose Hessian */
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(hessian);
    const unsigned int nEV = es.eigenvalues().size();
    if (modeToMaximize < 0) {
      throw std::runtime_error("The requested mode index " + std::to_string(modeToMaximize) +
                               " is negative, which is not possible.");
    }
    if (static_cast<unsigned>(modeToMaximize) >= nEV) {
      throw std::runtime_error("The requested mode index " + std::to_string(modeToMaximize) + " is too large for " +
                               std::to_string(nEV) + " eigenvectors.");
    }

    Eigen::VectorXd minstep = es.eigenvectors().transpose() * gradients;
    /* determine followed EV */
    Eigen::VectorXd overlaps = Eigen::VectorXd::Zero(nEV);
    if (_prevFollowedEigenvector.size() != 0) {
      /* estimate overlap of EVs to previous EV, which was followed, by cosine similarity */
      for (unsigned i = 0; i < nEV; ++i) {
        /* only evaluate EVs with negative eigenvalues */
        if (es.eigenvalues()[i] >= 0.0) {
          if (i == 0) {
            throw NoNegativeEigenValueException();
          }
          break;
        }
        overlaps.row(i) = (_prevFollowedEigenvector.transpose() * es.eigenvectors().col(i)) /
                          (_prevFollowedEigenvector.norm() * es.eigenvectors().col(i).norm());
      }
      /* use absolute values of cosine */
      Eigen::VectorXd absOverlaps = overlaps.cwiseAbs();
      absOverlaps.maxCoeff(&modeToMaximize);
    }
    _prevFollowedEigenvector = es.eigenvectors().col(modeToMaximize);
    /* Generate lambda P */
    Eigen::Matrix2d tmp1;
    tmp1 << es.eigenvalues()[modeToMaximize], minstep[modeToMaximize], minstep[modeToMaximize], 0;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2;
    es2.compute(tmp1);
    const double lambdaP = es2.eigenvalues()[1]; // corresponds to highest eigenvalue
    /* set of eigenvalues and steps without followed EV */
    Eigen::VectorXd notFollowingMinstep(nEV - 1);
    Eigen::VectorXd notFollowingEigenvalues(nEV - 1);
    bool skipped = false;
    for (int i = 0; i < static_cast<int>(nEV); ++i) {
      if (i != modeToMaximize) {
        /* followed mode is not inserted, therefore shift by one once the mode was skipped */
        notFollowingMinstep.row(skipped ? i - 1 : i) = minstep.row(i);
        notFollowingEigenvalues.row(skipped ? i - 1 : i) = es.eigenvalues().row(i);
      }
      else {
        skipped = true;
      }
    }
    /* Generate lambda N */
    Eigen::MatrixXd tmp2(nEV, nEV);
    tmp2.block(0, 0, nEV - 1, nEV - 1) = notFollowingEigenvalues.asDiagonal();
    tmp2.block(nEV - 1, 0, 1, nEV - 1) = notFollowingMinstep.transpose();
    tmp2.block(0, nEV - 1, nEV - 1, 1) = notFollowingMinstep;
    tmp2(nEV - 1, nEV - 1) = 0;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es3;
    es3.compute(tmp2);
    const double lambdaN = es3.eigenvalues()[0]; // corresponds to lowest eigenvalue
    /* Contribution from followed eigenvalue to step */
    Eigen::VectorXd steps =
        (-minstep[modeToMaximize] / (es.eigenvalues()[modeToMaximize] - lambdaP)) * es.eigenvectors().col(modeToMaximize);
    if (steps != steps) {
      throw std::runtime_error("Step size in optimizer is nan.");
    }
    /* Contribution from other eigenvalues to step */
    for (int i = 0; i < static_cast<int>(nEV); ++i) {
      if (i != modeToMaximize) {
        steps -= (minstep[i] / (es.eigenvalues()[i] - lambdaN)) * es.eigenvectors().col(i);
      }
    }
    if (steps.hasNaN()) {
      throw std::runtime_error("Step size in optimizer is nan.");
    }
    return steps;
  };

  /*! @brief Screening vector for constraining.
   *
   * An empty mask indicates all parameters may change. True values indicate
   * a parameter is not constrained. Must be either empty or have the same
   * length as the parameter vector.
   */
  Eigen::Matrix<bool, Eigen::Dynamic, 1> mask;

  /*! @brief saved iteration values to check for oscillations
   * @note setting to < 3 corresponds to turning the check off
   */
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

  /**
   * @brief Alter parameters by steps, zeroing out changes to constrained parameters
   *
   * @param parameters Function parameters
   * @param steps Eigen template expression of delta to parameters
   */
  template<typename Derived>
  void constrainedAdd(Eigen::VectorXd& parameters, const Eigen::DenseBase<Derived>& steps) {
    if (mask.size() == 0) {
      parameters += steps;
    }
    else {
      parameters += mask.select(steps, 0);
    }
  }

  Eigen::VectorXd constrainGradient(const Eigen::VectorXd& gradient) {
    if (mask.size() == 0) {
      return gradient;
    }

    return mask.select(gradient, 0);
  }

  //! The cycle number the optimize function starts counting with
  int _startCycle = 1;
  //! The current cycle number of the optimize function
  int _cycle;
  //! Last iteration function values for oscillation testing
  std::deque<double> _valueMemory;

 private:
  //! A vector of all registered observers.
  std::vector<std::function<void(const int&, const double&, const Eigen::VectorXd&)>> _observers;
  //! The previously followed eigenvector
  Eigen::VectorXd _prevFollowedEigenvector;
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
