/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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
class LBFGS : public Optimizer {
 public:
  static constexpr const char* lbfgsMaxm = "lbfgs_maxm";
  static constexpr const char* lbfgsLinesearch = "lbfgs_linesearch";
  static constexpr const char* lbfgsC1 = "lbfgs_c1";
  static constexpr const char* lbfgsC2 = "lbfgs_c2";
  static constexpr const char* lbfgsStepLength = "lbfgs_step_length";
  static constexpr const char* lbfgsUseTrustRadius = "lbfgs_use_trust_radius";
  static constexpr const char* lbfgsTrustRadius = "lbfgs_trust_radius";

  /// @brief Default constructor.
  LBFGS() = default;
  /**
   * @brief The main routine of the optimizer that carries out the actual optimization.
   *
   * @tparam UpdateFunction A lambda function with a void return value, and the arguments:\n
   *                        1. const Eigen::VectorXd& parameters\n
   *                        2. double& value\n
   *                        3. Eigen::VectorXd& gradients
   *
   * @param parameters The parameters to be optimized.
   * @param function   The function to be evaluated in order to get values and gradients
   *                   for a given set of parameters.
   * @param check      The ConvergenceCheck to be used in order to determine when the optimization
   *                   is finished or should stop for other reasons.
   * @return int       Returns the number of optimization cycles carried out until the conclusion
   *                   of the optimization function.
   */
  template<class UpdateFunction>
  int optimize(Eigen::VectorXd& parameters, UpdateFunction&& function, GradientBasedCheck& check) {
    /* number of parameters treated */
    unsigned int nParams = parameters.size();
    double value;
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
    /* start with one steepest descent step */
    Eigen::VectorXd stepVector = -0.1 * stepLength * gradients;
    bool stop = false;
    int cycle = 0;
    while (!stop) {
      cycle++;
      /* Rotate the (now) old parameter values to the local variables. */
      parametersOld = parameters;
      /* Rotate the (now) old gradient values to the local variable. */
      gradientsOld = gradients;
      /* Rotate the (now) old value to the local variable. */
      oldValue = value;
      /* Check trust radius */
      if (useTrustRadius) {
        double rms = sqrt((stepVector * stepLength).squaredNorm() / stepVector.size());
        if (rms > trustRadius) {
          stepVector *= trustRadius / rms;
        }
      }
      /* Update the parameters */
      parameters.noalias() += (stepLength * stepVector);
      /* Update gradients/value and check convergence */
      function(parameters, value, gradients);
      this->triggerObservers(cycle, value, parameters);
      stop = check.checkMaxIterations(cycle);
      if (!stop) {
        stop = check.checkConvergence(parameters, value, gradients);
      }
      if (stop)
        break;

      /*===========================*
       *  Generate new stepVector
       *===========================*/
      /* Armijo condition (Wolfe condition I),
       * used keep the step length sufficient
       */
      bool armijo = value <= (oldValue + c1 * stepLength * gradientsOld.dot(stepVector));
      /* Curvature condition (Wolfe condition II) */
      bool curvature = gradients.dot(stepVector) >= (c2 * gradientsOld.dot(stepVector));
      /* Backtracking condition */
      bool backtracking = value > oldValue;

      /* Linesearch */
      if (linesearch) {
        if (!armijo) {
          /* Check need for backtracking */
          if (backtracking) {
            parameters = parametersOld;
            gradients = gradientsOld;
            /* Adjust step length */
            stepLength *= 0.5;
            continue;
          }
          /* Technically the step size is wrong and there should be backtracking
           * followed by a step with the correct length either way.
           * However, if the step was sane (i.e. the value fell), the step
           * is kept in order to save computation time on updates.
           */
          /* Adjust step length */
          stepLength *= 0.5;
        }
        else if (armijo && !curvature) {
          /* Check need for backtracking */
          if (backtracking) {
            parameters = parametersOld;
            gradients = gradientsOld;
            /* Adjust step length */
            stepLength *= 1.5;
            continue;
          }
          /* Technically the step size is wrong and there should be backtracking
           * followed by a step with the correct length either way.
           * However, if the step was sane (i.e. the value fell), the step
           * is kept in order to save computation time on updates.
           */
          /* Adjust step length */
          stepLength *= 1.5;
        }
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
          alpha[i] = dx.col(i).dot(stepVector) / 1.0e-6;
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
  };
  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include the L-BFGS's options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  virtual void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
    UniversalSettings::IntDescriptor lbfgs_maxm("The maximum number of old steps stored.");
    lbfgs_maxm.setDefaultValue(maxm);
    collection.push_back(LBFGS::lbfgsMaxm, lbfgs_maxm);
    UniversalSettings::BoolDescriptor lbfgs_linesearch("Switch to turn on and off the use of a linesearch.");
    lbfgs_linesearch.setDefaultValue(linesearch);
    collection.push_back(LBFGS::lbfgsLinesearch, lbfgs_linesearch);
    UniversalSettings::DoubleDescriptor lbfgs_c1("1st parameter for the Wolfe conditions.");
    lbfgs_c1.setDefaultValue(c1);
    collection.push_back(LBFGS::lbfgsC1, lbfgs_c1);
    UniversalSettings::DoubleDescriptor lbfgs_c2("2nd parameter for the Wolfe conditions.");
    lbfgs_c2.setDefaultValue(c2);
    collection.push_back(LBFGS::lbfgsC2, lbfgs_c2);
    UniversalSettings::DoubleDescriptor lbfgs_stepLength("The initial step length used in the L-BFGS.");
    lbfgs_stepLength.setDefaultValue(stepLength);
    collection.push_back(LBFGS::lbfgsStepLength, lbfgs_stepLength);
    UniversalSettings::BoolDescriptor lbfgs_useTrustRadius("Enable the use of a trust radius.");
    lbfgs_useTrustRadius.setDefaultValue(useTrustRadius);
    collection.push_back(LBFGS::lbfgsUseTrustRadius, lbfgs_useTrustRadius);
    UniversalSettings::DoubleDescriptor lbfgs_trustRadius("The maximum size (RMS) of a taken step.");
    lbfgs_trustRadius.setDefaultValue(trustRadius);
    collection.push_back(LBFGS::lbfgsTrustRadius, lbfgs_trustRadius);
  };

  /**
   * @brief Updates the L-BFGS's options with those values given in the Settings.
   * @param settings The settings to update the option of the steepest descent with.
   */
  virtual void applySettings(const Settings& settings) final {
    maxm = (unsigned int)(settings.getInt(LBFGS::lbfgsMaxm));
    linesearch = settings.getBool(LBFGS::lbfgsLinesearch);
    c1 = settings.getDouble(LBFGS::lbfgsC1);
    c2 = settings.getDouble(LBFGS::lbfgsC2);
    stepLength = settings.getDouble(LBFGS::lbfgsStepLength);
    useTrustRadius = settings.getBool(LBFGS::lbfgsUseTrustRadius);
    trustRadius = settings.getDouble(LBFGS::lbfgsTrustRadius);
  };
  /**
   * @brief The maximum number of old steps stored.
   *
   * After hitting this maximum number of stored data,
   * the oldest parameters and gradients will be replaced
   */
  unsigned int maxm = 50;
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
  /// @brief Enable the use of a trust radius.
  bool useTrustRadius = false;
  /// @brief The maximum size (RMS) of a taken step.
  double trustRadius = 0.1;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_LBFGS_H_
