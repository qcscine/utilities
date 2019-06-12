/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_BOFILL_H_
#define UTILS_BOFILL_H_

#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/Optimizer.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
namespace Scine {
namespace Utils {

/**
 * @brief An implementation of the Bofill optimization algorithm for saddlepoints.
 *
 * Minimizes along all coordinates except one, along which it maximizes.
 *
 * Implemented, as described in: Phys. Chem. Chem. Phys., 2002, 4, 11–15
 * The original paper by Bofill is the following:  J. Comput. Chem., 1994, 15, 1
 */
class Bofill : public Optimizer {
 public:
  static constexpr const char* bofillTrustRadius = "bofill_trust_radius";

  /// @brief Default constructor.
  Bofill() = default;
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
    /* Setting all gradients to zero. */
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(nParams);
    /* Initial update */
    function(parameters, value, gradients);
    this->triggerObservers(0, value, parameters);
    /* Initialize 'old' variables */
    Eigen::VectorXd xOld(parameters);
    Eigen::VectorXd gOld(gradients);
    /* Initialize Hessian matrix */
    Eigen::MatrixXd B = Eigen::MatrixXd::Identity(nParams, nParams);
    if (eigenvector) {
      parameters -= 0.01 * (gradients - 2.0 * (gradients.dot(*eigenvector) / eigenvector->norm()) * (*eigenvector));
    }
    else {
      parameters -= 0.01 * gradients;
    }
    function(parameters, value, gradients);
    bool stop = false;
    int cycle = 0;
    while (!stop) {
      cycle++;
      Eigen::VectorXd dx = parameters - xOld;
      Eigen::VectorXd dg = gradients - gOld;
      // Update Hessian
      const double tmp1 = dx.transpose() * dx;
      Eigen::VectorXd tmp2(dg - B * dx);
      const double tmpdotdx = tmp2.dot(dx);
      //  The Bofill weight factor
      double bofillFacor = (tmpdotdx * tmpdotdx) / (tmp2.dot(tmp2) * dx.dot(dx));
      //  Powell  symmetric Broyden (PSB) part of the Bofill algorithm
      B += (1 - bofillFacor) * ((tmp2 * dx.transpose() + dx * tmp2.transpose()) / tmp1);
      B -= (1 - bofillFacor) * ((tmpdotdx * dx * dx.transpose()) / (tmp1 * tmp1));
      //  SR1 part of the Bofill algorithm
      B += bofillFacor * ((tmp2 * tmp2.transpose()) / tmpdotdx);
      // Hessian decomposition
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B);
      Eigen::VectorXd minstep = (es.eigenvectors().transpose() * gradients).eval();
      // Generate lambda P
      Eigen::Matrix2d tmp3;
      tmp3 << es.eigenvalues()[0], minstep[0], minstep[0], 0;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2;
      es2.compute(tmp3);
      const double lambdaP = es2.eigenvalues()[1];
      // Generate lambda N
      const unsigned int nEV = es.eigenvalues().size();
      Eigen::MatrixXd tmp4 = Eigen::MatrixXd::Zero(nEV, nEV);
      tmp4.block(0, 0, nEV - 1, nEV - 1) = es.eigenvalues().tail(nEV - 1).asDiagonal();
      tmp4.block(nEV - 1, 0, 1, nEV - 1) = minstep.tail(nEV - 1).transpose();
      tmp4.block(0, nEV - 1, nEV - 1, 1) = minstep.tail(nEV - 1);
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es3;
      es3.compute(tmp4);
      const double lambdaN = es3.eigenvalues()[0];
      // Contribution from first eigenvalue to step
      Eigen::VectorXd steps = (-minstep(0) / (es.eigenvalues()[0] - lambdaP)) * es.eigenvectors().col(0);
      // Contribution from other eigenvalues to step
      for (int i = 1; i < es.eigenvalues().size(); ++i) {
        steps -= (minstep(i) / (es.eigenvalues()[i] - lambdaN)) * es.eigenvectors().col(i);
      }
      // // Check trustradius and take a step
      xOld = parameters;
      double rms = sqrt(steps.squaredNorm() / steps.size());
      if (rms > trustRadius) {
        steps *= trustRadius / rms;
      }
      parameters += steps;
      // Update
      gOld = gradients;
      function(parameters, value, gradients);
      // Check convergence
      this->triggerObservers(cycle, value, parameters);
      stop = check.checkMaxIterations(cycle);
      if (!stop) {
        stop = check.checkConvergence(parameters, value, gradients);
      }
    }
    return cycle;
  };
  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include the Bofill's options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  virtual void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
    UniversalSettings::DoubleDescriptor bofill_trust_radius("The maximum RMS of a taken step.");
    bofill_trust_radius.setDefaultValue(trustRadius);
    collection.push_back(Bofill::bofillTrustRadius, bofill_trust_radius);
  };
  /**
   * @brief Updates the Bofill's options with those values given in the Settings.
   * @param settings The settings to update the option of the steepest descent with.
   */
  virtual void applySettings(const Settings& settings) final {
    trustRadius = settings.getDouble(Bofill::bofillTrustRadius);
  };
  /// @brief The maximum RMS of a taken step.
  double trustRadius = 0.1;
  /// @brief The (optional) eigenvector for the first step.
  std::unique_ptr<Eigen::VectorXd> eigenvector;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_BOFILL_H_
