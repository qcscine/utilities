/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_LBFGS_H_
#define UTILS_LBFGS_H_

#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/Optimizer.h"
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
  static constexpr const char* lbfgsMaxm = "lbfgs_maxm";
  static constexpr const char* lbfgsLinesearch = "lbfgs_linesearch";
  static constexpr const char* lbfgsC1 = "lbfgs_c1";
  static constexpr const char* lbfgsC2 = "lbfgs_c2";
  static constexpr const char* lbfgsStepLength = "lbfgs_step_length";
  static constexpr const char* lbfgsUseTrustRadius = "lbfgs_use_trust_radius";
  static constexpr const char* lbfgsTrustRadius = "lbfgs_trust_radius";
  static constexpr const char* lbfgsMaxBacktracking = "lbfgs_max_backtracking";

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
  int optimize(Eigen::VectorXd& parameters, UpdateFunction&& function, ConvergenceCheckClass& check) {
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
    int cycle = _startCycle;
    this->triggerObservers(cycle, value, parameters);
    // Init parameters and value stored in the convergence checker
    check.setParametersAndValue(parameters, value);
    /* start with one steepest descent step */
    Eigen::VectorXd stepVector = -0.1 * stepLength * gradients;
    double currentStepLength = stepLength;
    bool stop = false;
    int btc = 0; // a counter for the number of consecutive backtracking steps
    while (!stop) {
      cycle++;
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
      this->triggerObservers(cycle, value, parameters);
      stop = check.checkMaxIterations(cycle) || check.checkConvergence(parameters, value, constrainGradient(gradients));
      // Check oscillation, perform new calculation if oscillating
      if (!stop && this->isOscillating(value)) {
        stepVector -= 0.5 * currentStepLength * stepVector;
        parametersOld.noalias() = parameters;
        gradientsOld.noalias() = gradients;
        oscillationCorrection(currentStepLength * stepVector, parameters);
        function(parameters, value, gradients);
        cycle++;
        this->triggerObservers(cycle, value, parameters);
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
    return cycle;
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
    collection.push_back(Lbfgs::lbfgsMaxm, lbfgs_maxm);
    UniversalSettings::BoolDescriptor lbfgs_linesearch("Switch to turn on and off the use of a linesearch.");
    lbfgs_linesearch.setDefaultValue(linesearch);
    collection.push_back(Lbfgs::lbfgsLinesearch, lbfgs_linesearch);
    UniversalSettings::DoubleDescriptor lbfgs_c1("1st parameter for the Wolfe conditions.");
    lbfgs_c1.setDefaultValue(c1);
    collection.push_back(Lbfgs::lbfgsC1, lbfgs_c1);
    UniversalSettings::DoubleDescriptor lbfgs_c2("2nd parameter for the Wolfe conditions.");
    lbfgs_c2.setDefaultValue(c2);
    collection.push_back(Lbfgs::lbfgsC2, lbfgs_c2);
    UniversalSettings::DoubleDescriptor lbfgs_stepLength("The initial step length used in the L-BFGS.");
    lbfgs_stepLength.setMinimum(0.0);
    lbfgs_stepLength.setDefaultValue(stepLength);
    collection.push_back(Lbfgs::lbfgsStepLength, lbfgs_stepLength);
    UniversalSettings::BoolDescriptor lbfgs_useTrustRadius("Enable the use of a trust radius for all steps.");
    lbfgs_useTrustRadius.setDefaultValue(useTrustRadius);
    collection.push_back(Lbfgs::lbfgsUseTrustRadius, lbfgs_useTrustRadius);
    UniversalSettings::DoubleDescriptor lbfgs_trustRadius("The maximum size (RMS) of a taken step.");
    lbfgs_trustRadius.setMinimum(0.0);
    lbfgs_trustRadius.setDefaultValue(trustRadius);
    collection.push_back(Lbfgs::lbfgsTrustRadius, lbfgs_trustRadius);
    UniversalSettings::IntDescriptor lbfgs_maxBacktracking(
        "The maximum number of consecutive backtracking steps allowed.");
    lbfgs_maxBacktracking.setMinimum(0);
    lbfgs_maxBacktracking.setDefaultValue(maxBacktracking);
    collection.push_back(Lbfgs::lbfgsMaxBacktracking, lbfgs_maxBacktracking);
  };

  /**
   * @brief Updates the L-BFGS's options with those values given in the Settings.
   * @param settings The settings to update the option of the L-BFGS with.
   */
  void applySettings(const Settings& settings) final {
    maxm = (unsigned int)(settings.getInt(Lbfgs::lbfgsMaxm));
    linesearch = settings.getBool(Lbfgs::lbfgsLinesearch);
    c1 = settings.getDouble(Lbfgs::lbfgsC1);
    c2 = settings.getDouble(Lbfgs::lbfgsC2);
    stepLength = settings.getDouble(Lbfgs::lbfgsStepLength);
    useTrustRadius = settings.getBool(Lbfgs::lbfgsUseTrustRadius);
    trustRadius = settings.getDouble(Lbfgs::lbfgsTrustRadius);
    maxBacktracking = settings.getInt(Lbfgs::lbfgsMaxBacktracking);
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
