/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_EIGENVECTORFOLLOWING_H_
#define UTILS_EIGENVECTORFOLLOWING_H_

#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/Optimizer.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>

namespace Scine {
namespace Utils {

/**
 * @brief An implementation of an Eigenvector following optimization algorithm.
 *
 * This algorithm is intended to find the maximum along one single eigenvector
 * and the minimum along all other eigenvectors of a given system/Hessian.
 */
class EigenvectorFollowing : public Optimizer {
 public:
  static constexpr const char* evTrustRadius = "ev_trust_radius";
  /// @brief Default constructor.
  EigenvectorFollowing() = default;
  /**
   * @brief The main routine of the optimizer that carries out the actual optimization.
   *
   * @tparam UpdateFunction A lambda function with a void return value, and the arguments:\n
   *                        1. const Eigen::VectorXd& parameters\n
   *                        2. double& value\n
   *                        3. Eigen::VectorXd& gradients
   *                        4. Eigen::Matrixd& the Hessian
   *                        5. bool a flag if the Hessian is to be calculated
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
    /* Initialize gradients and Hessian. */
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(nParams);
    auto hessian = std::make_unique<Eigen::MatrixXd>(nParams, nParams);
    hessian->setZero();
    /* Initial update */
    function(parameters, value, gradients, *hessian, true);
    this->triggerObservers(0, value, parameters);
    bool stop = false;
    int cycle = 0;
    Eigen::VectorXd prevFollowedEigenvector;
    prevFollowedEigenvector.resize(nParams);
    while (!stop) {
      cycle++;
      /* Decompose Hessian */
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
      es.compute(*hessian);
      const unsigned int nEV = es.eigenvalues().size();
      Eigen::VectorXd minstep = es.eigenvectors().transpose() * gradients;
      /* determine followed EV */
      /* for first step lowest eigenvalue is followed; could be changed to user-defined value */
      unsigned int modeToFollow = 0;
      Eigen::VectorXd overlaps = Eigen::VectorXd::Zero(nEV);
      if (cycle != 1) {
        /* estimate overlap of EVs to previous EV, which was followed, by cosine similarity */
        for (int i = 0; i < nEV; ++i) {
          /* only evaluate EVs with negative eigenvalues */
          if (es.eigenvalues()[i] < 0) {
            overlaps.row(i) = (prevFollowedEigenvector.transpose() * es.eigenvectors().col(i)) /
                              (prevFollowedEigenvector.norm() * es.eigenvectors().col(i).norm());
          }
          else
            break;
        }
        /* use absolute values of cosine */
        Eigen::VectorXd absOverlaps = overlaps.cwiseAbs();
        if (absOverlaps.size() > 0) {
          absOverlaps.maxCoeff(&modeToFollow);
        }
        else {
          std::cout
              << "WARNING: Did not find negative eigenvalue, trying to maximize along mode of the lowest eigenvalue."
              << std::endl;
          modeToFollow = 0;
        }
      }
      prevFollowedEigenvector = es.eigenvectors().col(modeToFollow);
      /* Generate lambda P */
      Eigen::Matrix2d tmp1;
      tmp1 << es.eigenvalues()[modeToFollow], minstep[modeToFollow], minstep[modeToFollow], 0;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2;
      es2.compute(tmp1);
      const double lambdaP = es2.eigenvalues()[1]; // corresponds to highest eigenvalue
      /* set of eigenvalues and steps without followed EV */
      Eigen::VectorXd notFollowingMinstep(nEV - 1);
      Eigen::VectorXd notFollowingEigenvalues(nEV - 1);
      bool skipped = false;
      for (int i = 0; i < nEV; ++i) {
        if (i != modeToFollow) {
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
          (-minstep[modeToFollow] / (es.eigenvalues()[modeToFollow] - lambdaP)) * es.eigenvectors().col(modeToFollow);
      /* Contribution from other eigenvalues to step */
      for (int i = 0; i < nEV; ++i) {
        if (i != modeToFollow) {
          steps -= (minstep[i] / (es.eigenvalues()[i] - lambdaN)) * es.eigenvectors().col(i);
        }
      }
      /* Take a step */
      double rms = sqrt(steps.squaredNorm() / steps.size());
      if (rms > trustRadius) {
        steps *= trustRadius / rms;
      }
      parameters += steps;
      // Update and check convergence
      function(parameters, value, gradients, *hessian, true);
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
   *        thus expanding it to include the optimizers's options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  virtual void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
    UniversalSettings::DoubleDescriptor ev_trust_radius("The maximum RMS of a taken step.");
    ev_trust_radius.setDefaultValue(trustRadius);
    collection.push_back(EigenvectorFollowing::evTrustRadius, ev_trust_radius);
  };
  /**
   * @brief Updates theoptimizers's options with those values given in the Settings.
   * @param settings The settings to update the option of the steepest descent with.
   */
  virtual void applySettings(const Settings& settings) final {
    trustRadius = settings.getDouble(EigenvectorFollowing::evTrustRadius);
  };
  /// @brief The maximum RMS of a taken step.
  double trustRadius = 0.3;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_EIGENVECTORFOLLOWING_H_
