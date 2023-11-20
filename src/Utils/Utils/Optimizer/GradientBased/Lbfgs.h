/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_LBFGS_H_
#define UTILS_LBFGS_H_

#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/Optimizer.h"
#include "Utils/UniversalSettings/OptimizationSettingsNames.h"
#include <Core/Log.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @brief An implementation of a limited memory BFGS optimization algorithm.
 *
 * The implementation includes an optional linesearch using both Wolfe conditions.
 */
class Lbfgs : public Optimizer {
 public:
  /// @brief Default constructor.
  Lbfgs() = default;
  /**
   * @brief The main routine of the optimizer that carries out the actual optimization.
   *
   * @tparam UpdateFunction A lambda function with a void return value, and the arguments:\n
   *                        1. const Eigen::VectorXd& parameters\n
   *                        2. double& value\n
   *                        3. Eigen::VectorXd& gradients
   * @tparam ConvergenceCheckClass The convergence check class, needs to have a check function signature
   *                               that is the same as the GradientBasedCheck.
   *
   * @param parameters The parameters to be optimized.
   * @param function   The function to be evaluated in order to get values and gradients
   *                   for a given set of parameters.
   * @param check      The ConvergenceCheck to be used in order to determine when the optimization
   *                   is finished or should stop for other reasons.
   * @return int       Returns the number of optimization cycles carried out until the conclusion
   *                   of the optimization function.
   */
  template<class UpdateFunction, class ConvergenceCheckClass>
  int optimize(Eigen::VectorXd& parameters, UpdateFunction&& function, ConvergenceCheckClass& check,
               Core::Log& /* log */) {
    /* number of parameters treated */
    unsigned int nParams = parameters.size();
    if (nParams == 0) {
      throw EmptyOptimizerParametersException();
    }
    double value = 0.0;
    unsigned int m = 0;

    Eigen::MatrixXd dg(nParams, maxm);
    dg.setZero();
    Eigen::MatrixXd dx(nParams, maxm);
    dx.setZero();

    /* Setting all gradients to zero. */
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(nParams);
    /* Setting all old gradients to zero. */
    Eigen::VectorXd gradientsOld = Eigen::VectorXd::Zero(nParams);
    /* Set "assumed previous parameters" as parametersOld. */
    Eigen::VectorXd parametersOld(parameters);
    double oldValue = value;
    function(parameters, value, gradients);
    _cycle = _startCycle;
    this->triggerObservers(_cycle, value, parameters);
    // Init parameters and value stored in the convergence checker
    check.setParametersAndValue(parameters, value);
    /* start with one steepest descent step */
    Eigen::VectorXd stepVector = -0.1 * stepLength * gradients;
    double currentStepLength = stepLength;
    bool stop = false;
    int btc = 0; // a counter for the number of consecutive backtracking steps
    while (!stop) {
      _cycle++;
      /* Rotate the (now) old parameter values to the local variables. */
      parametersOld.noalias() = parameters;
      /* Rotate the (now) old gradient values to the local variable. */
      gradientsOld.noalias() = gradients;
      /* Rotate the (now) old value to the local variable. */
      oldValue = value;
      /* Check trust radius */
      if (useTrustRadius) {
        double maxVal = stepVector.array().abs().maxCoeff();
        if (maxVal * currentStepLength > trustRadius) {
          stepVector *= trustRadius / (maxVal * currentStepLength);
        }
      }
      /* Update the parameters */
      constrainedAdd(parameters, currentStepLength * stepVector);
      /* Update gradients/value and check convergence */
      function(parameters, value, gradients);
      this->triggerObservers(_cycle, value, parameters);
      stop = check.checkMaxIterations(_cycle) || check.checkConvergence(parameters, value, constrainGradient(gradients));
      // Check oscillation, perform new calculation if oscillating
      if (!stop && this->isOscillating(value)) {
        stepVector -= 0.5 * currentStepLength * stepVector;
        parametersOld.noalias() = parameters;
        gradientsOld.noalias() = gradients;
        oscillationCorrection(currentStepLength * stepVector, parameters);
        function(parameters, value, gradients);
        _cycle++;
        this->triggerObservers(_cycle, value, parameters);
      }
      if (stop) {
        break;
      }

      /*===========================*
       *  Generate new stepVector
       *===========================*/
      /* Armijo condition (Wolfe condition I),
       * used keep the step length sufficient
       */
      bool armijo = value <= (oldValue + c1 * currentStepLength * gradientsOld.dot(stepVector));
      /* Curvature condition (Wolfe condition II) */
      bool curvature = gradients.dot(stepVector) >= (c2 * gradientsOld.dot(stepVector));

      /* Linesearch */
      if (linesearch) {
        if ((btc > maxBacktracking) && (!armijo || (armijo && !curvature))) {
          // If there was backtracking for more than 5 consecutive steps,
          // and there would be another backtracking step now, then:
          //  backtrack, but clear all stored data, and start with a
          //  fresh gradient descent
          parameters = parametersOld;
          gradients = gradientsOld;
          value = oldValue;
          currentStepLength = stepLength;
          stepVector = -0.1 * stepLength * gradients;
          btc = 0;
          m = 0;
          dg.setZero();
          dx.setZero();
          continue;
        }
        if (!armijo) {
          /* Check need for backtracking */
          parameters = parametersOld;
          gradients = gradientsOld;
          value = oldValue;
          /* Adjust step length */
          currentStepLength *= 0.523;
          btc++;
          continue;
        }

        if (armijo && !curvature) {
          /* Check need for backtracking */
          parameters = parametersOld;
          gradients = gradientsOld;
          value = oldValue;
          /* Adjust step length */
          currentStepLength *= 1.515;
          btc++;
          continue;
        }

        btc = 0;
      }

      /* Normal L-BFGS update, with possibly adjusted step size. */
      stepVector = -gradients;
      if (m < maxm) {
        dg.col(m) = gradients - gradientsOld;
        dx.col(m) = parameters - parametersOld;
        ++m;
      }
      else {
        dg.leftCols(maxm - 1) = dg.rightCols(maxm - 1);
        dx.leftCols(maxm - 1) = dx.rightCols(maxm - 1);
        dg.col(maxm - 1) = gradients - gradientsOld;
        dx.col(maxm - 1) = parameters - parametersOld;
      }
      /* The actual L-BFGS update */
      Eigen::VectorXd alpha(m);
      for (int i = m - 1; i > -1; --i) {
        double dxDotdg = dx.col(i).dot(dg.col(i));
        if (fabs(dxDotdg) < 1.0e-6) {
          alpha[i] = dx.col(i).dot(stepVector) / ((dxDotdg < 0.0) ? -1.0e-6 : 1.0e-6);
        }
        else {
          alpha[i] = dx.col(i).dot(stepVector) / dxDotdg;
        }
        stepVector.noalias() -= alpha[i] * dg.col(i);
      }
      stepVector *= (dx.col(m - 1).dot(dg.col(m - 1)) / dg.col(m - 1).squaredNorm());
      for (unsigned int i = 0; i < m; ++i) {
        double dxDotdg = dx.col(i).dot(dg.col(i));
        double beta = dg.col(i).dot(stepVector);
        if (fabs(dxDotdg) < 1.0e-6) {
          beta /= 1.0e-6;
        }
        else {
          beta /= dx.col(i).dot(dg.col(i));
        }
        stepVector.noalias() += (alpha[i] - beta) * dx.col(i);
      }
    }
    return _cycle;
  }
  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include the L-BFGS's options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
    UniversalSettings::IntDescriptor lbfgs_maxm("The maximum number of old steps stored.");
    lbfgs_maxm.setMinimum(0);
    lbfgs_maxm.setDefaultValue(maxm);
    collection.push_back(SettingsNames::Optimizations::Lbfgs::maxm, lbfgs_maxm);
    UniversalSettings::BoolDescriptor _linesearch("Switch to turn on and off the use of a linesearch.");
    _linesearch.setDefaultValue(linesearch);
    collection.push_back(SettingsNames::Optimizations::Lbfgs::linesearch, _linesearch);
    UniversalSettings::DoubleDescriptor _c1("1st parameter for the Wolfe conditions.");
    _c1.setDefaultValue(c1);
    collection.push_back(SettingsNames::Optimizations::Lbfgs::c1, _c1);
    UniversalSettings::DoubleDescriptor _c2("2nd parameter for the Wolfe conditions.");
    _c2.setDefaultValue(c2);
    collection.push_back(SettingsNames::Optimizations::Lbfgs::c2, _c2);
    UniversalSettings::DoubleDescriptor _stepLength("The initial step length used in the L-BFGS.");
    _stepLength.setMinimum(0.0);
    _stepLength.setDefaultValue(stepLength);
    collection.push_back(SettingsNames::Optimizations::Lbfgs::stepLength, _stepLength);
    UniversalSettings::BoolDescriptor _useTrustRadius("Enable the use of a trust radius for all steps.");
    _useTrustRadius.setDefaultValue(useTrustRadius);
    collection.push_back(SettingsNames::Optimizations::Lbfgs::useTrustRadius, _useTrustRadius);
    UniversalSettings::DoubleDescriptor _trustRadius("The maximum size (RMS) of a taken step.");
    _trustRadius.setMinimum(0.0);
    _trustRadius.setDefaultValue(trustRadius);
    collection.push_back(SettingsNames::Optimizations::Lbfgs::trustRadius, _trustRadius);
    UniversalSettings::IntDescriptor _maxBacktracking("The maximum number of consecutive backtracking steps allowed.");
    _maxBacktracking.setMinimum(0);
    _maxBacktracking.setDefaultValue(maxBacktracking);
    collection.push_back(SettingsNames::Optimizations::Lbfgs::maxBacktracking, _maxBacktracking);
  };

  /**
   * @brief Updates the L-BFGS's options with those values given in the Settings.
   * @param settings The settings to update the option of the L-BFGS with.
   */
  void applySettings(const Settings& settings) final {
    maxm = (unsigned int)(settings.getInt(SettingsNames::Optimizations::Lbfgs::maxm));
    linesearch = settings.getBool(SettingsNames::Optimizations::Lbfgs::linesearch);
    c1 = settings.getDouble(SettingsNames::Optimizations::Lbfgs::c1);
    c2 = settings.getDouble(SettingsNames::Optimizations::Lbfgs::c2);
    stepLength = settings.getDouble(SettingsNames::Optimizations::Lbfgs::stepLength);
    useTrustRadius = settings.getBool(SettingsNames::Optimizations::Lbfgs::useTrustRadius);
    trustRadius = settings.getDouble(SettingsNames::Optimizations::Lbfgs::trustRadius);
    maxBacktracking = settings.getInt(SettingsNames::Optimizations::Lbfgs::maxBacktracking);
    if (!useTrustRadius && std::fabs(trustRadius - 0.1) > 1e-6) {
      throw std::logic_error("A trust radius was specified, but the trust radius was not activated. "
                             "Please also set the setting 'lbfgs_use_trust_radius': true, if you specify a radius.");
    }
  };
  /**
   * @brief The maximum number of old steps stored.
   *
   * After hitting this maximum number of stored data,
   * the oldest parameters and gradients will be replaced
   */
  unsigned int maxm = 10;
  /// @brief Switch to turn on and off the use of a linesearch
  bool linesearch = true;
  /**
   * @brief 1st parameter for the Wolfe conditions.
   *
   * This parameter is only relevant if the linesearch is turned on.
   */
  double c1 = 0.0001;
  /**
   * @brief 2nd parameter for the Wolfe conditions.
   *
   * This parameter is only relevant if the linesearch is turned on.
   * Also a value as low as 0.1 can be used.
   */
  double c2 = 0.9;
  /**
   * @brief The initial step length used in the L-BFGS.
   *
   * Note: the first step is a gradient descent with 0.1 times the steplength.
   */
  double stepLength = 1.0;
  /// @brief Enable the use of a trust radius for all steps.
  bool useTrustRadius = false;
  /// @brief The maximum size (RMS) of a taken step.
  double trustRadius = 0.1;
  /// @brief The maximum number of consecutive backtracking steps allowed.
  int maxBacktracking = 5;
  /// @brief TODO A possible Hessian projection
  // std::function<void(Eigen::MatrixXd&)> projection;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_LBFGS_H_
