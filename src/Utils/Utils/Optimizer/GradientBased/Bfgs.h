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
#include "Utils/UniversalSettings/OptimizationSettingsNames.h"
#include <Core/Log.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @brief An implementation of the BFGS optimization algorithm.
 *
 * The implementation includes an optional GDIIS.
 *
 * The initial approximate inverse Hessian matrix in the BFGS (cycle > 1) is build by multiplying the identity matrix I
 * by dxTdg / dgTdg. This is an attempt to make the size of the first approximate inverse Hessian similar to the real
 * inverse Hessian. See in:
 * (2006) Quasi-Newton Methods. In: Numerical Optimization by J. Nocedal, S. J. Wright
 * [DOI: 10.1007/978-0-387-40065-5_6], https://link.springer.com/chapter/10.1007/978-0-387-40065-5_6
 *
 * The implementation includes a damping of the BFGS update. The damping should keep the inverse approximate Hessian
 * positive definite. The idea origins from Powell (1978), but the implementation is inspired by a slightly modified
 * version (in the reference called BFGS1, sigma2 is 0.9, sigma3 is 9.0):
 * Adv. Model. Optim., 2009, 11, 1, 63
 * https://camo.ici.ro/journal/vol11/v11a5.pdf
 *
 * As the algorithm works with the approximate inverse Hessian, the update is performed in the inverse space as well.
 * The connection is reported in:
 * arXiv:2006.08877 [cs.LG]
 * https://arxiv.org/abs/2006.08877
 *
 */
class Bfgs : public Optimizer {
 public:
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
  int optimize(Eigen::VectorXd& parameters, UpdateFunction&& function, ConvergenceCheckClass& check,
               Core::Log& /* log */) {
    /* number of parameters treated */
    unsigned int nParams = parameters.size();
    if (nParams == 0) {
      throw EmptyOptimizerParametersException();
    }
    if (check.maxIter <= static_cast<unsigned int>(minIter)) {
      throw std::logic_error(
          "Maximum number of cycles is smaller or equal to the minimum number of cycles.\n"
          "Please set the maximum number of cycles at least one larger than the minimum number of cycles.");
    }
    double value = 0.0;

    /* Initialize Hessian inverse if needed */
    if (invH.size() == 0) {
      invH = 0.5 * Eigen::MatrixXd::Identity(nParams, nParams);
    }
    /* Initialize GDIIS */
    Gdiis diis(invH, gdiisMaxStore);
    /* Setting all gradients to zero. */
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(nParams);
    /* Setting all old gradients to zero. */
    Eigen::VectorXd gradientsOld = Eigen::VectorXd::Zero(nParams);
    /* Set "assumed previous parameters" as parametersOld. */
    Eigen::VectorXd parametersOld(parameters);
    /* Get start cycle */
    _cycle = _startCycle;
    function(parameters, value, gradients);
    this->triggerObservers(_cycle, value, parameters);
    // Init parameters and value stored in the convergence checker
    check.setParametersAndValue(parameters, value);

    /* start with one steepest descent step */
    Eigen::VectorXd stepVector = -invH * gradients;
    /* Rotate the (now) old parameter values to the local variables. */
    parametersOld.noalias() = parameters;
    /* Rotate the (now) old gradient values to the local variable. */
    gradientsOld.noalias() = gradients;
    /* Check trust radius for first SD Step */
    if (useTrustRadius) {
      Eigen::VectorXd projectedStep = parameters + stepVector - parametersOld;
      double maxVal = projectedStep.array().abs().maxCoeff();
      if (maxVal > trustRadius) {
        /* Scale step vector to trust radius */
        stepVector.noalias() = projectedStep * trustRadius / maxVal;
      }
    }
    /* Update the parameters */
    constrainedAdd(parameters, stepVector);
    bool stop = false;
    bool minIterPerformed = false;
    while (!stop) {
      _cycle++;
      function(parameters, value, gradients);
      this->triggerObservers(_cycle, value, parameters);
      /* Once minIterPerformed is true, it no longer needs to be checked */
      if (!minIterPerformed) {
        minIterPerformed = (_cycle > minIter) ? true : false;
      }
      /* Check stop criteria only after requested initial updates have been performed */
      stop = (check.checkMaxIterations(_cycle) || check.checkConvergence(parameters, value, constrainGradient(gradients))) &&
             minIterPerformed;
      if (stop) {
        break;
      }

      /* BFGS inverse Hessian update */
      Eigen::VectorXd dx = parameters - parametersOld;
      Eigen::VectorXd dg = gradients - gradientsOld;
      double dxTdg = dx.dot(dg);
      const double dgTinvHdg = dg.transpose() * invH * dg;
      /* Set initial inverse Hessian to dxTdg / dgTdg * I */
      if (_cycle == _startCycle + 1) {
        invH.diagonal() *= 2.0 * dxTdg / dg.dot(dg);
      }

      /* Powell Update from Al-Baali Gandetti */
      const double sigma2 = 0.9;
      const double sigma3 = 9.0;
      double delta = 1.0;
      if (fabs(dxTdg) < fabs((1.0 - sigma2) * dgTinvHdg)) {
        delta = sigma2 * dgTinvHdg / (dgTinvHdg - dxTdg);
      }
      else if (fabs(dxTdg) > fabs((1.0 + sigma3) * dgTinvHdg)) {
        delta = -sigma3 * dgTinvHdg / (dgTinvHdg - dxTdg);
      }
      /* Update dx by dx = delta * dx + (1.0 - delta) * invH * dg */
      if (fabs(delta - 1.0) > 1e-12) {
        dx.noalias() = delta * dx + (1.0 - delta) * invH * dg;
        dxTdg = dx.dot(dg);
      }
      /* Sanity check for dxTdg */
      if (fabs(dxTdg) < 1e-9) {
        dxTdg = (dxTdg < 0.0) ? -1e-9 : 1e-9;
      }

      if (dxTdg > 0.0) {
        /* Update inverse approximate Hessian */
        const double alpha = (dxTdg + dg.transpose() * invH * dg) / (dxTdg * dxTdg);
        const double beta = 1.0 / dxTdg;
        invH += alpha * (dx * dx.transpose()) - beta * (invH * dg * dx.transpose() + dx * dg.transpose() * invH);
      }
      /* Reset invH, if positive definitness is not given*/
      else {
        invH.noalias() = Eigen::MatrixXd::Identity(nParams, nParams) * 0.5;
        /* Flush GDIIS */
        if (useGdiis) {
          diis.flush();
        }
      }
      if (projection) {
        projection(invH);
      }

      /* Rotate the (now) old parameter values to the local variables. */
      parametersOld.noalias() = parameters;
      /* Rotate the (now) old gradient values to the local variable. */
      gradientsOld.noalias() = gradients;
      /* Update the parameters */
      if (useGdiis) {
        diis.update(parameters, gradients);
        if (mask.size() > 0) {
          /* Revert constrained parameters to old values */
          parameters = mask.select(parameters, parametersOld);
        }
      }
      stepVector.noalias() = -invH * gradients;

      /* Check trust radius */
      /* Steps will be scaled by trust radius divided by the maximum value in the projected step */
      if (useTrustRadius) {
        Eigen::VectorXd projectedStep = parameters + stepVector - parametersOld;
        double maxVal = projectedStep.array().abs().maxCoeff();
        if (maxVal > trustRadius) {
          /* Scale step vector to trust radius */
          stepVector.noalias() = projectedStep * trustRadius / maxVal;
          /* Reset parameters back to before gdiis */
          parameters.noalias() = parametersOld;
          /* Reset inverse Hessian matrix and the GDIIS */
          this->resetInverseHessian(dg, dxTdg);
          if (useGdiis) {
            diis.flush();
          }
        }
      }
      constrainedAdd(parameters, stepVector);
      if (this->isOscillating(value)) {
        this->oscillationCorrection(stepVector, parameters);
        this->resetInverseHessian(dg, dxTdg);
        this->_valueMemory.clear();
        if (useGdiis) {
          diis.flush();
        }
      }
    }
    return _cycle;
  }
  /**
   * @brief Reset the inverse hessian.
   * @param dg    Gradient change between iterations.
   * @param dxTdg Scalar product of parameter and gradient change.
   */
  void resetInverseHessian(const Eigen::VectorXd& dg, const double& dxTdg) {
    const double dgTdg = dg.transpose() * dg;
    const unsigned int nParams = dg.rows();
    if (dgTdg > 1e-9) {
      invH.noalias() = Eigen::MatrixXd::Identity(nParams, nParams) * dxTdg / dgTdg;
    }
    else {
      invH.noalias() = Eigen::MatrixXd::Identity(nParams, nParams) * 0.5;
    }
  }

  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include the BFGS's options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
    UniversalSettings::IntDescriptor bfgs_minIter(
        "The minimal number of cycles to be performed before the stop criteria is checked.");
    bfgs_minIter.setDefaultValue(minIter);
    bfgs_minIter.setMinimum(1.0);
    collection.push_back(SettingsNames::Optimizations::Bfgs::minIter, bfgs_minIter);
    UniversalSettings::BoolDescriptor bfgs_useTrustRadius("Enable the use of a trust radius for all steps.");
    bfgs_useTrustRadius.setDefaultValue(useTrustRadius);
    collection.push_back(SettingsNames::Optimizations::Bfgs::useTrustRadius, bfgs_useTrustRadius);
    UniversalSettings::DoubleDescriptor bfgs_trustRadius("The maximum size (RMS) of a taken step.");
    bfgs_trustRadius.setMinimum(0.0);
    bfgs_trustRadius.setDefaultValue(trustRadius);
    collection.push_back(SettingsNames::Optimizations::Bfgs::trustRadius, bfgs_trustRadius);
    UniversalSettings::BoolDescriptor bfgs_useGdiis(
        "Switch to enable the use of a GDIIS possibly accelerating convergence");
    bfgs_useGdiis.setDefaultValue(useGdiis);
    collection.push_back(SettingsNames::Optimizations::Bfgs::useGdiis, bfgs_useGdiis);
    UniversalSettings::IntDescriptor bfgs_gdiisMaxStore("The maximum number of old steps used in the GDIIS.");
    bfgs_gdiisMaxStore.setMinimum(0);
    bfgs_gdiisMaxStore.setDefaultValue(gdiisMaxStore);
    collection.push_back(SettingsNames::Optimizations::Bfgs::gdiisMaxStore, bfgs_gdiisMaxStore);
  };

  /**
   * @brief Updates the BFGS's options with those values given in the Settings.
   * @param settings The settings to update the option of the BFGS with.
   */
  void applySettings(const Settings& settings) final {
    minIter = settings.getInt(SettingsNames::Optimizations::Bfgs::minIter);
    useTrustRadius = settings.getBool(SettingsNames::Optimizations::Bfgs::useTrustRadius);
    trustRadius = settings.getDouble(SettingsNames::Optimizations::Bfgs::trustRadius);
    useGdiis = settings.getBool(SettingsNames::Optimizations::Bfgs::useGdiis);
    gdiisMaxStore = settings.getInt(SettingsNames::Optimizations::Bfgs::gdiisMaxStore);
    if (!useTrustRadius && std::fabs(trustRadius - 0.3) > 1e-6) {
      throw std::logic_error("A trust radius was specified, but the trust radius was not activated. "
                             "Please also set the setting 'bfgs_use_trust_radius': true, if you specify a radius.");
    }
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
  void prepareRestart(const int cycleNumber) final {
    /* Set start cycle number */
    _startCycle = cycleNumber;
    /* Remove inverse Hessian and projection function */
    (this->invH).resize(0, 0);
    this->projection = nullptr;
    /* Clear value memory for oscillating correction */
    _valueMemory.clear();
  };
  /**
   * @brief Minimal number of cycles before the stop criteria are checked.
   *
   * The B matrix gets initialised in cycle 1,
   * which corresponds to an SD step with the gradients scaled by 0.5.
   * In cycle 2, the diagonal elements get set to certain values and the matrix updated once.
   *
   * Per default, with minIter set to 1, in cycle 2 the stop criteria are checked before anything
   * is done with the matrix.
   * Hence, the optimization can converge with 0 updates on the B matrix .
   *
   * Example for the case minIter is set to 2:
   * The stop criteria are not checked in cycle 2 and the B matrix gets manipulated as described above.
   * In cycle 3, the stop criteria are checked again.
   */
  int minIter = 1;
  /// @brief Enable the use of a trust radius for all steps.
  bool useTrustRadius = true;
  /// @brief The maximum size (RMS) of a taken step.
  double trustRadius = 0.3;
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
