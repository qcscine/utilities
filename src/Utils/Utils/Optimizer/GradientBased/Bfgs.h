/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_BFGS_H_
#define UTILS_BFGS_H_

#include "Utils/Optimizer/GradientBased/Gdiis.h"
#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/Optimizer.h"
#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @brief An implementation of the BFGS optimization algorithm.
 *
 * The implementation includes an optional GDIIS.
 */
class Bfgs : public Optimizer {
 public:
  static constexpr const char* bfgsUseTrustRadius = "bfgs_use_trust_radius";
  static constexpr const char* bfgsTrustRadius = "bfgs_trust_radius";
  static constexpr const char* bfgsUseGdiis = "bfgs_use_gdiis";
  static constexpr const char* bfgsGdiisMaxStore = "bfgs_gdiis_max_store";

  /// @brief Default constructor.
  Bfgs() = default;
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

    /* Initialize Hessian inverse if needed */
    if (invH.size() == 0) {
      invH = 0.1 * Eigen::MatrixXd::Identity(nParams, nParams);
    }
    /* Initialize GDIIS */
    Gdiis diis(invH, gdiisMaxStore);
    /* Setting all gradients to zero. */
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(nParams);
    /* Setting all old gradients to zero. */
    Eigen::VectorXd gradientsOld = Eigen::VectorXd::Zero(nParams);
    /* Set "assumed previous parameters" as parametersOld. */
    Eigen::VectorXd parametersOld(parameters);
    double oldValue = value;
    /* Get start cycle */
    int cycle = _startCycle;
    function(parameters, value, gradients);
    this->triggerObservers(cycle, value, parameters);
    // Init parameters and value stored in the convergence checker
    check.setParametersAndValue(parameters, value);
    /* start with one steepest descent step */
    Eigen::VectorXd stepVector = -invH * gradients;
    bool stop = false;
    while (!stop) {
      cycle++;
      /* Rotate the (now) old parameter values to the local variables. */
      parametersOld.noalias() = parameters;
      /* Rotate the (now) old gradient values to the local variable. */
      gradientsOld.noalias() = gradients;
      /* Rotate the (now) old value to the local variable. */
      oldValue = value;

      if ((sqrt(gradients.squaredNorm() / gradients.size()) < 0.2) && useGdiis && cycle - _startCycle > gdiisMaxStore) {
        parameters = diis.update(parameters, gradients);
        stepVector.noalias() = (parameters - parametersOld);
        /* Check trust radius */
        if (useTrustRadius) {
          double rms = sqrt(stepVector.squaredNorm() / stepVector.size());
          if (rms > trustRadius) {
            parameters.noalias() += ((trustRadius / rms) - 1.0) * stepVector;
          }
        }
        function(parameters, value, gradients);
        this->triggerObservers(cycle, value, parameters);
      }
      else {
        /* Check trust radius */
        if (useTrustRadius) {
          double rms = sqrt((stepVector).squaredNorm() / stepVector.size());
          if (rms > trustRadius) {
            stepVector *= trustRadius / rms;
          }
        }
        /* Update the parameters */
        parameters.noalias() += stepVector;
        /* Update gradients/value and check convergence */
        function(parameters, value, gradients);
        this->triggerObservers(cycle, value, parameters);
        if (useGdiis)
          diis.store(parameters, gradients);
      }
      if (value > oldValue) {
        /* Backtrack of value rises */
        parameters.noalias() = parametersOld;
        gradients.noalias() = gradientsOld;
        invH = Eigen::MatrixXd::Identity(nParams, nParams);
        if (useGdiis)
          diis.flush();
      }

      stop = check.checkMaxIterations(cycle);
      if (!stop) {
        stop = check.checkConvergence(parameters, value, gradients);
      }
      if (stop)
        break;
      // Check oscillation, perform new calculation if oscillating
      if (this->isOscillating(value)) {
        parametersOld.noalias() = parameters;
        gradientsOld.noalias() = gradients;
        this->oscillationCorrection(stepVector, parameters);
        function(parameters, value, gradients);
        cycle++;
        this->triggerObservers(cycle, value, parameters);
      }
      /* BFGS inverse Hessian update */
      Eigen::VectorXd dx = parameters - parametersOld;
      Eigen::VectorXd dg = gradients - gradientsOld;
      const double dxTdx = dx.dot(dx);
      double dxTdg = dx.dot(dg);
      if (fabs(dxTdg) < 1e-9)
        dxTdg = (dxTdg < 0.0) ? -1e-9 : 1e-9;
      const double alpha = (dxTdg + dg.transpose() * invH * dg) / (dxTdg * dxTdg);
      const double beta = 1.0 / dxTdg;
      invH += alpha * (dx * dx.transpose()) - beta * (invH * dg * dx.transpose() + dx * dg.transpose() * invH);
      if (projection) {
        projection(invH);
      }
      stepVector.noalias() = -invH * gradients;
    }
    return cycle;
  };
  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include the BFGS's options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  virtual void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
    UniversalSettings::BoolDescriptor bfgs_useTrustRadius("Enable the use of a trust radius for all steps.");
    bfgs_useTrustRadius.setDefaultValue(useTrustRadius);
    collection.push_back(Bfgs::bfgsUseTrustRadius, bfgs_useTrustRadius);
    UniversalSettings::DoubleDescriptor bfgs_trustRadius("The maximum size (RMS) of a taken step.");
    bfgs_trustRadius.setDefaultValue(trustRadius);
    collection.push_back(Bfgs::bfgsTrustRadius, bfgs_trustRadius);
    UniversalSettings::BoolDescriptor bfgs_useGdiis(
        "Switch to enable the use of a GDIIS possibly accelerating convergence");
    bfgs_useGdiis.setDefaultValue(useGdiis);
    collection.push_back(Bfgs::bfgsUseGdiis, bfgs_useGdiis);
    UniversalSettings::IntDescriptor bfgs_gdiisMaxStore("The maximum number of old steps used in the GDIIS.");
    bfgs_gdiisMaxStore.setDefaultValue(gdiisMaxStore);
    collection.push_back(Bfgs::bfgsGdiisMaxStore, bfgs_gdiisMaxStore);
  };

  /**
   * @brief Updates the BFGS's options with those values given in the Settings.
   * @param settings The settings to update the option of the BFGS with.
   */
  virtual void applySettings(const Settings& settings) final {
    useTrustRadius = settings.getBool(Bfgs::bfgsUseTrustRadius);
    trustRadius = settings.getDouble(Bfgs::bfgsTrustRadius);
    useGdiis = settings.getBool(Bfgs::bfgsUseGdiis);
    gdiisMaxStore = settings.getInt(Bfgs::bfgsGdiisMaxStore);
  };
  /**
   * @brief Prepares the BFGS optimizer for rerunning its optimize function.
   *
   * This function is used to prepare an BFGS optimizer instance for rerunning
   * its optimize function with possibly different settings.
   * It changes the cycle count the optimizer starts with when the optimize
   * function is called and removes the stored inverse Hessian and its
   * projection function.
   * The value memory for an eventual oscillating correction is cleared.
   *
   * @param cycleNumber The cycle number the optimizer starts with.
   */
  virtual void prepareRestart(const int& cycleNumber) final {
    // Set start cycle number
    _startCycle = cycleNumber;
    // Remove inverse Hessian and projection function
    (this->invH).resize(0, 0);
    this->projection = nullptr;
    // Clear value memory for oscillating correction
    _initializedValueMemory = false;
    _valueMemory.clear();
  };
  /// @brief Enable the use of a trust radius for all steps.
  bool useTrustRadius = false;
  /// @brief The maximum size (RMS) of a taken step.
  double trustRadius = 0.1;
  /// @brief Switch to enable the use of a GDIIS possibly accelerating convergence.
  bool useGdiis = true;
  /// @brief The maximum number of old steps used in the GDIIS.
  int gdiisMaxStore = 5;
  /// @brief The inverse Hessian
  Eigen::MatrixXd invH;
  /// @brief A possible Hessian projection
  std::function<void(Eigen::MatrixXd&)> projection;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_BFGS_H_
