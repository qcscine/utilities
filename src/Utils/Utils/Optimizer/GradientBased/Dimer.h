/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_DIMER_H_
#define UTILS_DIMER_H_

#include "Utils/Geometry/InternalCoordinates.h"
#include "Utils/Optimizer/GradientBased/Gdiis.h"
#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/Optimizer.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
namespace Scine {
namespace Utils {

/**
 * @brief An implementation of the Dimer optimization algorithm for saddlepoints.
 *
 * transition state search without hessian by finding minimum curvature by gradient calculations of a dimer on PES
 * dimer is constructed by estimated eigenvector or random vector, then dimer is rotated iteratively until it is aligned
 * with lowest curvature mode, then translation is performed
 * this is repeated until converged
 *
 * Implemented, as described in: Kaestner J. and Sherwood P., Journal of Chemical Physics, 2008
 *                               Shang C. and Liu ZP., Journal of Chemical Theory and Computation, 2010
 */
class Dimer : public Optimizer {
 public:
  static constexpr const char* dimerSkipFirstRotation = "dimer_skip_first_rotation";
  static constexpr const char* dimerDecreaseRotationGradientThreshold = "dimer_decrease_rotation_gradient_threshold";
  static constexpr const char* dimerGradientInterpolation = "dimer_gradient_interpolation";
  static constexpr const char* dimerRotationCG = "dimer_rotation_conjugate_gradient";
  static constexpr const char* dimerRotationLBFGS = "dimer_rotation_lbfgs";
  static constexpr const char* dimerOnlyOneRotation = "dimer_only_one_rotation";
  static constexpr const char* dimerUseGdiis = "dimer_gdiis";
  static constexpr const char* dimerUseProjectionTrustRadius = "dimer_projection_trust_radius";
  static constexpr const char* dimerMultiScale = "dimer_multi_scale";
  static constexpr const char* dimerRadius = "dimer_radius";
  static constexpr const char* dimerPhiTolerance = "dimer_phi_tolerance";
  static constexpr const char* dimerRotationGradientThresholdFirstCycle = "dimer_rotation_gradient_first";
  static constexpr const char* dimerRotationGradientThresholdOtherCycles = "dimer_rotation_gradient_other";
  static constexpr const char* dimerLoweredRotationGradientThreshold = "dimer_lowered_rotation_gradient";
  static constexpr const char* dimerGradRMSDthreshold = "dimer_grad_rmsd_threshold";
  static constexpr const char* dimerTrustRadius = "dimer_trust_radius";
  static constexpr const char* dimerDefaultTranslationStep = "dimer_default_translation_step";
  static constexpr const char* dimerMaxRotationsFirstCycle = "dimer_max_rotations_first_cycle";
  static constexpr const char* dimerMaxRotationsOtherCycles = "dimer_max_rotations_other_cycle";
  static constexpr const char* dimerIntervalOfRotations = "dimer_interval_of_rotations";
  static constexpr const char* dimerCycleOfRotationGradientDecrease = "dimer_cycle_of_rotation_gradient_decrease";
  static constexpr const char* dimerMaxBacktracking = "dimer_max_backtracking";
  /// @brief Default constructor.
  Dimer() = default;
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
    setVectorsToZero(nParams);
    /* Construct GDIIS object with identity matrix */
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(nParams, nParams);
    Gdiis gdiis(identity, maxBacktracking);
    bool previousMaxScaling = false;
    bool appliedStepsizeScaling = false;
    double value, value1, curvature;
    unsigned int m = 0;
    double currentStepLength = defaultTranslationStep;
    sanityCheck();
    createDimerAxis(nParams);
    /* First calculation */
    function(parameters, value, _gradients);
    this->triggerObservers(0, value, parameters);
    unsigned int cycle = 0;
    bool stop = false;
    while (!stop) {
      cycle++;
      /* Generate Dimer */
      _parametersR1.noalias() = parameters + radius * _dimerAxis;
      /* Force parallel to the dimer axis */
      _fParaOld.noalias() = _fPara;
      _fPara.noalias() = -_gradients.dot(_dimerAxis) * _dimerAxis;
      /* check whether rotation is performed in this cycle */
      int maxRotationCycles = (cycle == 1) ? maxRotationsFirstCycle : maxRotationsOtherCycles;
      if (determineIfPerformRotation(cycle, parameters, function)) {
        _numberOfPerformedRotationCycles++;
        curvature = (cycle == 1) ? rotationWithPhi(parameters, value, maxRotationCycles, function)
                                 : rotationWithGradient(parameters, value, maxRotationCycles, function);
        /* Reset scaling after rotation */
        currentStepLength = defaultTranslationStep;
      }
      if (curvature > 0) {
        /* no negative eigenvalue estimate --> still in convex region, trying to solely maximize lowest mode */
        _modForce.noalias() = _gradients.dot(_dimerAxis) * _dimerAxis;
      }
      else
        _modForce.noalias() = 2.0 * (_gradients.dot(_dimerAxis) * _dimerAxis) - _gradients;
      /* G corresponds to step vector */
      _G.noalias() = _modForce;
      Eigen::VectorXd negG = -_modForce;
      /* if modification to step size, gradient is checked and scaling is activated below certain threshold */
      if (useGdiis || useProjectionTrustRadius) {
        /* scaling has been activated, no need to calculate gradientRMSD, otherwise calculate */
        if (appliedStepsizeScaling || (std::sqrt(_gradients.squaredNorm() / _gradients.size())) < gradRMSDthreshold) {
          appliedStepsizeScaling = true; // bool turned true if entered by gradient calculations
          /* steps determined by gdiis object */
          if (useGdiis) {
            parameters = gdiis.update(parameters, negG); // directly returns parameters
          }
          /* projection method does not change direction and length has been changed at the end of last cycle */
          else if (useProjectionTrustRadius) {
            _steps = currentStepLength * _G;
          }
        }
        /* scaling cannot be applied yet, but gdiis can store data and both methods perform steepest descent */
        else {
          if (useGdiis) {
            gdiis.store(parameters, negG);
          }
          _steps = defaultTranslationStep * _G;
        }
      }
      /* no modification, steps determined by unchanged steplength and vector */
      else {
        _steps = defaultTranslationStep * _G;
      }
      if (!useGdiis || !appliedStepsizeScaling) {
        /* Check trustradius and take a step */
        double rms = sqrt(_steps.squaredNorm() / _steps.size());
        if (rms > trustRadius) {
          _steps *= trustRadius / rms;
        }
        parameters += _steps;
      }
      /* Calculate energy and gradient of new position */
      function(parameters, value, _gradients);
      this->triggerObservers(cycle, value, parameters);
      /* Check convergence */
      stop = check.checkMaxIterations(cycle);
      if (!stop)
        stop = check.checkConvergence(parameters, value, _gradients);
      if (useProjectionTrustRadius && appliedStepsizeScaling && !stop)
        currentStepLength = scaleWithProjectionTrustRadius(parameters, value, function, previousMaxScaling);
    } // end of optimization while loop
    std::cout << "\n    Number of individual rotations: " << _sumRotations << std::endl;
    return cycle;
  };

  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include the Dimers's options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  virtual void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
    UniversalSettings::BoolDescriptor skip_first_rotation(
        "If the provided guessVector shall be used without rotation in first step.");
    skip_first_rotation.setDefaultValue(skipFirstRotation);
    collection.push_back(Dimer::dimerSkipFirstRotation, skip_first_rotation);
    UniversalSettings::BoolDescriptor decrease_rotation_gradient_threshold(
        "If threshold for convergence in rotation shall be lowered after certain number of rotationcycles.");
    decrease_rotation_gradient_threshold.setDefaultValue(decreaseRotationGradientThreshold);
    collection.push_back(Dimer::dimerDecreaseRotationGradientThreshold, decrease_rotation_gradient_threshold);
    UniversalSettings::BoolDescriptor gradient_interpolation(
        "If gradient of R1 shall be estimated by interpolation from last rotation data.");
    gradient_interpolation.setDefaultValue(gradientInterpolation);
    collection.push_back(Dimer::dimerGradientInterpolation, gradient_interpolation);
    UniversalSettings::BoolDescriptor rotation_conjugate_gradient("If conjugate gradient shall be used for rotation.");
    rotation_conjugate_gradient.setDefaultValue(rotationCG);
    collection.push_back(Dimer::dimerRotationCG, rotation_conjugate_gradient);
    UniversalSettings::BoolDescriptor rotation_lbfgs("If L-BFGS shall be used for rotation.");
    rotation_lbfgs.setDefaultValue(rotationLBFGS);
    collection.push_back(Dimer::dimerRotationLBFGS, rotation_lbfgs);
    UniversalSettings::BoolDescriptor only_one_rotation("If rotation shall only be performed in first step.");
    only_one_rotation.setDefaultValue(onlyOneRotation);
    collection.push_back(Dimer::dimerOnlyOneRotation, only_one_rotation);
    UniversalSettings::BoolDescriptor use_gdiis("If GDIIS shall be used for translation.");
    use_gdiis.setDefaultValue(useGdiis);
    collection.push_back(Dimer::dimerUseGdiis, use_gdiis);
    UniversalSettings::BoolDescriptor use_projection_trust_radius(
        "If the step size shall be estimated by projection of the gradient at current step onto the last dimer axis.");
    use_projection_trust_radius.setDefaultValue(useProjectionTrustRadius);
    collection.push_back(Dimer::dimerUseProjectionTrustRadius, use_projection_trust_radius);
    UniversalSettings::BoolDescriptor multi_scale(
        "If the projection step size factor shall be applied to previous step size.");
    multi_scale.setDefaultValue(multiScale);
    collection.push_back(Dimer::dimerMultiScale, multi_scale);
    UniversalSettings::DoubleDescriptor dimer_radius("Distance between the two images of the dimer.");
    dimer_radius.setDefaultValue(radius);
    collection.push_back(Dimer::dimerRadius, dimer_radius);
    UniversalSettings::DoubleDescriptor phi_tolerance("The convergence criterion for the rotation angle.");
    phi_tolerance.setDefaultValue(phiTolerance);
    collection.push_back(Dimer::dimerPhiTolerance, phi_tolerance);
    UniversalSettings::DoubleDescriptor rotation_gradient_threshold_first_cycle(
        "The convergence criterion for the gradient of the rotation in the first rotation cycle.");
    rotation_gradient_threshold_first_cycle.setDefaultValue(rotationGradientThresholdFirstCycle);
    collection.push_back(Dimer::dimerRotationGradientThresholdFirstCycle, rotation_gradient_threshold_first_cycle);
    UniversalSettings::DoubleDescriptor rotation_gradient_threshold_other_cycles(
        "The convergence criterion for the gradient of the rotation in all but the first rotation cycle.");
    rotation_gradient_threshold_other_cycles.setDefaultValue(rotationGradientThresholdOtherCycles);
    collection.push_back(Dimer::dimerRotationGradientThresholdOtherCycles, rotation_gradient_threshold_other_cycles);
    UniversalSettings::DoubleDescriptor lowered_rotation_gradient_threshold(
        "The convergence criterion for the gradient of the rotation if the criterion is lowered after some cycles.");
    lowered_rotation_gradient_threshold.setDefaultValue(loweredRotationGradientThreshold);
    collection.push_back(Dimer::dimerLoweredRotationGradientThreshold, lowered_rotation_gradient_threshold);
    UniversalSettings::DoubleDescriptor grad_rmsd_threshold(
        "RMSD threshold of the gradient for starting the stepsize scaling.");
    grad_rmsd_threshold.setDefaultValue(gradRMSDthreshold);
    collection.push_back(Dimer::dimerGradRMSDthreshold, grad_rmsd_threshold);
    UniversalSettings::DoubleDescriptor dimer_trust_radius("The maximum RMS of a taken step.");
    dimer_trust_radius.setDefaultValue(trustRadius);
    collection.push_back(Dimer::dimerTrustRadius, dimer_trust_radius);
    UniversalSettings::DoubleDescriptor dimer_default_translation_step(
        "The factor to multiple the stepsize vector in a steepest descent step.");
    dimer_default_translation_step.setDefaultValue(defaultTranslationStep);
    collection.push_back(Dimer::dimerDefaultTranslationStep, dimer_default_translation_step);
    UniversalSettings::IntDescriptor max_rotations_first_cycle(
        "The maximum number of rotations in the first rotation cycle.");
    max_rotations_first_cycle.setDefaultValue(maxRotationsFirstCycle);
    collection.push_back(Dimer::dimerMaxRotationsFirstCycle, max_rotations_first_cycle);
    UniversalSettings::IntDescriptor max_rotations_other_cycle(
        "The maximum number of rotations in all but the first rotation cycle.");
    max_rotations_other_cycle.setDefaultValue(maxRotationsOtherCycles);
    collection.push_back(Dimer::dimerMaxRotationsOtherCycles, max_rotations_other_cycle);
    UniversalSettings::IntDescriptor interval_of_rotations(
        "The interval of performed rotation cycles in the total optimization steps.");
    interval_of_rotations.setDefaultValue(intervalOfRotations);
    collection.push_back(Dimer::dimerIntervalOfRotations, interval_of_rotations);
    UniversalSettings::IntDescriptor cycle_of_rotation_gradient_decrease(
        "The number of rotation cycles after which the threshold for the gradient of rotation is decreased.");
    cycle_of_rotation_gradient_decrease.setDefaultValue(cycleOfRotationGradientDecrease);
    collection.push_back(Dimer::dimerCycleOfRotationGradientDecrease, cycle_of_rotation_gradient_decrease);
    UniversalSettings::IntDescriptor max_back_tracking("The number of saved gradients during rotation in L-BFGS.");
    max_back_tracking.setDefaultValue(maxBacktracking);
    collection.push_back(Dimer::dimerMaxBacktracking, max_back_tracking);
  };
  /**
   * @brief Updates the Dimer's options with those values given in the Settings.
   * @param settings The settings to update the option of the steepest descent with.
   */
  virtual void applySettings(const Settings& settings) final {
    skipFirstRotation = settings.getBool(Dimer::dimerSkipFirstRotation);
    decreaseRotationGradientThreshold = settings.getBool(Dimer::dimerDecreaseRotationGradientThreshold);
    gradientInterpolation = settings.getBool(Dimer::dimerGradientInterpolation);
    rotationCG = settings.getBool(Dimer::dimerRotationCG);
    rotationLBFGS = settings.getBool(Dimer::dimerRotationLBFGS);
    onlyOneRotation = settings.getBool(Dimer::dimerOnlyOneRotation);
    useGdiis = settings.getBool(Dimer::dimerUseGdiis);
    useProjectionTrustRadius = settings.getBool(Dimer::dimerUseProjectionTrustRadius);
    multiScale = settings.getBool(Dimer::dimerMultiScale);
    radius = settings.getDouble(Dimer::dimerRadius);
    phiTolerance = settings.getDouble(Dimer::dimerPhiTolerance);
    rotationGradientThresholdFirstCycle = settings.getDouble(Dimer::dimerRotationGradientThresholdFirstCycle);
    rotationGradientThresholdOtherCycles = settings.getDouble(Dimer::dimerRotationGradientThresholdOtherCycles);
    loweredRotationGradientThreshold = settings.getDouble(Dimer::dimerLoweredRotationGradientThreshold);
    gradRMSDthreshold = settings.getDouble(Dimer::dimerGradRMSDthreshold);
    trustRadius = settings.getDouble(Dimer::dimerTrustRadius);
    defaultTranslationStep = settings.getDouble(Dimer::dimerDefaultTranslationStep);
    maxRotationsFirstCycle = settings.getInt(Dimer::dimerMaxRotationsFirstCycle);
    maxRotationsOtherCycles = settings.getInt(Dimer::dimerMaxRotationsOtherCycles);
    intervalOfRotations = settings.getInt(Dimer::dimerIntervalOfRotations);
    cycleOfRotationGradientDecrease = settings.getInt(Dimer::dimerCycleOfRotationGradientDecrease);
    maxBacktracking = settings.getInt(Dimer::dimerMaxBacktracking);
  };

  //// @brief If the provided guessVector shall be used without rotation in first step.
  bool skipFirstRotation = false;
  //// @brief If threshold for convergence in rotation shall be lowered after certain number of rotationcycles.
  bool decreaseRotationGradientThreshold = false;
  //// @brief If gradient of R1 shall be estimated by interpolation from last rotation data.
  bool gradientInterpolation = false;
  //// @brief If conjugate gradient shall be used for rotation.
  bool rotationCG = false;
  //// @brief If L-BFGS shall be used for rotation.
  bool rotationLBFGS = true;
  //// @brief If rotation shall only be performed in first step.
  bool onlyOneRotation = false;
  //// @brief If GDIIS shall be used for translation.
  bool useGdiis = false;
  //// @brief If the step size shall be estimated by projection of the gradient at current step onto the last dimer axis.
  bool useProjectionTrustRadius = true;
  //// @brief If the projection step size factor shall be applied to previous step size.
  bool multiScale = true;
  //// @brief Distance between the two images of the dimer
  double radius = 0.01;
  //// @brief The convergence criterion for the rotation angle.
  double phiTolerance = 1e-3;
  //// @brief The convergence criterion for the gradient of the rotation in the first rotation cycle.
  double rotationGradientThresholdFirstCycle = 1e-7;
  //// @brief The convergence criterion for the gradient of the rotation in all but the first rotation cycle.
  double rotationGradientThresholdOtherCycles = 1e-4;
  //// @brief The convergence criterion for the gradient of the rotation if the criterion is lowered after some cycles.
  double loweredRotationGradientThreshold = 1e-3;
  //// @brief The employed threshold depending on conditions.
  double rotationGradientThreshold = 0.0;
  //// @brief The convergence criterion for the gradient of the rotation if the criterion is lowered after some cycles.
  double gradRMSDthreshold = 1e-3;
  /// @brief The maximum RMS of a taken step.
  double trustRadius = 0.5;
  /// @brief The factor to multiple the stepsize vector in a steepest descent step.
  double defaultTranslationStep = 1;
  //// @brief The maximum number of rotations in the first rotation cycle.
  unsigned int maxRotationsFirstCycle = 100;
  //// @brief The maximum number of rotations in all but the first rotation cycle.
  unsigned int maxRotationsOtherCycles = 100;
  //// @brief The interval of performed rotation cycles in the total optimization steps.
  unsigned int intervalOfRotations = 5;
  //// @brief The number of rotation cycles after which the threshold for the gradient of rotation is decreased.
  unsigned int cycleOfRotationGradientDecrease = 5;
  //// @brief The number of saved gradients during rotation in L-BFGS.
  unsigned int maxBacktracking = 5;
  //// @brief The (optional) dimer axis for the first step provided by hessian calculation or trajectory.
  std::unique_ptr<Eigen::VectorXd> guessVector;

 private:
  void setVectorsToZero(const int& nParams) {
    _gradients = Eigen::VectorXd::Zero(nParams);
    _gradientsR1 = Eigen::VectorXd::Zero(nParams);
    _fPara = Eigen::VectorXd::Zero(nParams);
    _fParaOld = Eigen::VectorXd::Zero(nParams);
    _modForce = Eigen::VectorXd::Zero(nParams);
    _G = Eigen::VectorXd::Zero(nParams);
    _orthoFN = Eigen::VectorXd::Zero(nParams);
    _orthoG = Eigen::VectorXd::Zero(nParams);
    _previousOrthoFN = Eigen::VectorXd::Zero(nParams);
    _previousOrthoG = Eigen::VectorXd::Zero(nParams);
    _steps = Eigen::VectorXd::Zero(nParams);
    _dx.resize(nParams, maxBacktracking);
    _dg.resize(nParams, maxBacktracking);
    _dx.setZero();
    _dg.setZero();
  };

  void sanityCheck() {
    if (useGdiis && useProjectionTrustRadius) {
      /* Projection radius is default, so GDIIS is wished to be performed */
      useProjectionTrustRadius = false;
    }
    if (rotationLBFGS && rotationCG) {
      /* BFGS is default, so CG is wished to be performed */
      rotationLBFGS = false;
    }
  };

  /**
   * @brief Create first dimer
   *
   * If guessvector is given, it is used for the dimer axis
   * if no guess, a random vector is used
   *
   * @param nParams            The number of parameters to be optimized
   *
   * @return void
   */
  void createDimerAxis(const int& nParams) {
    bool randomNecessary = false;
    if (guessVector != nullptr) {
      if ((*guessVector).size() == nParams)
        _dimerAxis.noalias() = *guessVector;
      else {
        std::cout << "ERROR, given guess vector has different size than parameters. "
                     "Continuing transition state search, but omitting the guess."
                  << std::endl;
        guessVector = NULL;
        randomNecessary = true;
      }
    }
    else
      randomNecessary = true;
    if (randomNecessary) {
      srand(42);
      _dimerAxis.noalias() = Eigen::VectorXd::Random(nParams);
    }
    _dimerAxis.normalize();
  };

  /**
   * @brief Determine whether a rotation of the dimer shall be carried out in this optimization step
   *
   * @param cycle              The number of the current cycle
   * @param parameters         The parameters to be optimized
   * @tparam UpdateFunction    A lambda function with a void return value, and the arguments:\n
   *                           1. const Eigen::VectorXd& parameters\n
   *                           2. double& value\n
   *                           3. Eigen::VectorXd& gradients
   *
   * @return bool              Returns whether rotation shall be carried out
   */
  template<class UpdateFunction>
  bool determineIfPerformRotation(const int& cycle, const Eigen::VectorXd& parameters, UpdateFunction&& function) {
    double value1;
    double curvature;
    if (cycle == 1) {
      if (skipFirstRotation)
        return false;
      else {
        rotationGradientThreshold = rotationGradientThresholdFirstCycle;
        function(_parametersR1, value1, _gradientsR1);
        /* determine ortho force at current position to rotate */
        _orthoFN = 2.0 * ((_gradientsR1 - _gradients).dot(_dimerAxis)) * _dimerAxis - 2.0 * (_gradientsR1 - _gradients);
        return true;
      }
    }
    /* if only first rotation wants to be performed */
    else if (onlyOneRotation)
      return false;
    else if (cycle % intervalOfRotations == 0) {
      rotationGradientThreshold = rotationGradientThresholdOtherCycles;
      function(_parametersR1, value1, _gradientsR1);
      curvature = (_gradientsR1 - _gradients).dot(_dimerAxis) / radius;
      /* determine ortho force at current position to rotate */
      _orthoFN = 2.0 * ((_gradientsR1 - _gradients).dot(_dimerAxis)) * _dimerAxis - 2.0 * (_gradientsR1 - _gradients);
      if ((curvature > 0 && _fPara.norm() < _fParaOld.norm()) || (curvature < 0 && _fPara.norm() > _fParaOld.norm()))
        return true;
      else if (_orthoFN.norm() > 1e-3)
        return true;
    }
    return false;
  };

  /**
   * @brief Determine step vector of rotation with the Conjugate Gradient (Polak Ribiere) method
   *
   * @param rotationCycles     The number of individual rotation in the current rotation - 1
   *
   * @return Eigen::VectorXd   Returns step vector of rotation
   */
  Eigen::VectorXd conjugateGradientRotation(const int& rotationCycles) {
    _orthoG = _orthoFN;
    if (rotationCycles != 0) {
      /* Polak Ribiere */
      double beta = _orthoFN.dot(_orthoFN - _previousOrthoFN) / (_previousOrthoFN.dot(_previousOrthoFN));
      _orthoG = _orthoFN + beta * _previousOrthoG;
    }
    return _orthoG;
  };

  /**
   * @brief Determine the step vector of rotation based on the L-BFGS methods
   *
   * @param rotationCycles     The number of individual rotation in the current rotation - 1
   * @param m                  The counter for the L-BFGS backtracking
   *
   * @return Eigen::VectorXd   Returns step vector of rotation
   */
  Eigen::VectorXd lbfgsRotation(const int& rotationCycles, unsigned int& m) {
    _orthoG = _orthoFN;
    int maxm = maxBacktracking;
    if (rotationCycles == 0) {
      _dx.setZero();
      _dg.setZero();
      return _orthoG;
    }
    /* Normal L-BFGS update, with possibly adjusted step size. */
    if (m < maxm) {
      _dg.col(m) = _previousOrthoFN - _orthoFN;
      _dx.col(m) = _parametersR1 - _oldParametersR1;
      ++m;
    }
    else {
      _dg.leftCols(maxm - 1) = _dg.rightCols(maxm - 1);
      _dx.leftCols(maxm - 1) = _dx.rightCols(maxm - 1);
      _dg.col(maxm - 1) = _previousOrthoFN - _orthoFN;
      _dx.col(maxm - 1) = _parametersR1 - _oldParametersR1;
    }
    /* The actual L-BFGS update */
    Eigen::VectorXd alpha(m);
    for (int i = m - 1; i > -1; --i) {
      double dxDotdg = _dx.col(i).dot(_dg.col(i));
      if (fabs(dxDotdg) < 1.0e-6) {
        alpha[i] = _dx.col(i).dot(_orthoG) / ((dxDotdg < 0.0) ? -1.0e-6 : 1.0e-6);
      }
      else {
        alpha[i] = _dx.col(i).dot(_orthoG) / dxDotdg;
      }
      _orthoG.noalias() -= alpha[i] * _dg.col(i);
    }
    _orthoG *= (_dx.col(m - 1).dot(_dg.col(m - 1)) / _dg.col(m - 1).squaredNorm());
    for (unsigned int i = 0; i < m; ++i) {
      double dxDotdg = _dx.col(i).dot(_dg.col(i));
      double beta = _dg.col(i).dot(_orthoG);
      if (fabs(dxDotdg) < 1.0e-6) {
        beta /= 1.0e-6;
      }
      else {
        beta /= _dx.col(i).dot(_dg.col(i));
      }
      _orthoG.noalias() += (alpha[i] - beta) * _dx.col(i);
    }
    return _orthoG;
  };

  /**
   * @brief Flips dimeraxis if energy decreases along axis
   *
   * Rotate dimer axis if necessary, so dimer axis is afterwards correctly aligned with the direction in which the mode
   * can be optimized
   *
   * @param value              The value of the optimized point
   * @param parameters         The parameters to be optimized
   * @param value1             The value at R1
   * @param curvature          The curvature along the dimer axis
   * @tparam UpdateFunction    A lambda function with a void return value, and the arguments:\n
   *                           1. const Eigen::VectorXd& parameters\n
   *                           2. double& value\n
   *                           3. Eigen::VectorXd& gradients
   *
   * @return void
   */
  template<class UpdateFunction>
  void determineDirectionOfDimerToBeMaximized(const double& value, const Eigen::VectorXd& parameters, double& value1,
                                              double& curvature, UpdateFunction&& function) {
    if (value > value1) {
      _dimerAxis *= -1;
      _parametersR1.noalias() = parameters + _dimerAxis * radius;
      function(_parametersR1, value1, _gradientsR1);
      curvature = (_gradientsR1 - _gradients).dot(_dimerAxis) / radius;
    }
  };

  /**
   * @brief rotate dimer with trial angle
   *
   * Rotate the dimer based on a trial rotation in a determined plane and then determine angle of minimum curvature
   * requires 2 calculations per rotation and might not converge fast, but is able to find global minimum more reliable
   *
   * @param parameters         The parameters to be optimized
   * @param value              The value of the optimized point
   * @param maxRotationCycles  The number of maximum allowed rotations
   * @tparam UpdateFunction    A lambda function with a void return value, and the arguments:\n
   *                           1. const Eigen::VectorXd& parameters\n
   *                           2. double& value\n
   *                           3. Eigen::VectorXd& gradients
   *
   * @return double            Returns the curvature estimate of the rotated dimer
   */
  template<class UpdateFunction>
  double rotationWithPhi(const Eigen::VectorXd& parameters, const double& value, const int& maxRotationCycles,
                         UpdateFunction&& function) {
    bool reachedTolerance = false;
    unsigned int m = 0; // for lbfgs
    double cPhiMin = 0.0;
    double phiMin = 0.0;
    double phi1 = 0.0;
    double value1 = 0.0;
    Eigen::VectorXd gradientsR1phi1 = Eigen::VectorXd::Zero(parameters.size());
    /* starting individual rotations */
    for (int rotationCycles = 0; rotationCycles < maxRotationCycles; ++rotationCycles) {
      /* calculate E and F at R1 */
      /* either estimate gradient by interpolation or calculate it */
      if (gradientInterpolation && rotationCycles != 0) {
        _gradientsR1 = std::sin(phi1 - phiMin) / std::sin(phi1) * _gradientsR1 +
                       std::sin(phiMin) / std::sin(phi1) * gradientsR1phi1 +
                       (1 - std::cos(phiMin) - std::sin(phiMin) * std::tan(phi1 / 2)) * _gradients;
      }
      else {
        function(_parametersR1, value1, _gradientsR1);
      }
      /* determine ortho force at current position to determin rotation plane */
      Eigen::VectorXd orthogonalDimerAxis = _dimerAxis; // used for transformation of the dimer axis
      _orthoFN.noalias() =
          2.0 * ((_gradientsR1 - _gradients).dot(_dimerAxis)) * _dimerAxis - 2.0 * (_gradientsR1 - _gradients);
      _orthoG.noalias() = _orthoFN;
      /* Conjugate Gradient option for rotation */
      if (rotationCG)
        _orthoG = conjugateGradientRotation(rotationCycles);
      /* L-BFGS for rotation */
      else if (rotationLBFGS)
        _orthoG = lbfgsRotation(rotationCycles, m);
      _orthoG -= _orthoG.dot(_dimerAxis) * _dimerAxis; // ensures that modified direction is still orthogonal to dimer
                                                       // axis
      /* unit vector which spans rotation plane together with dimer axis */
      Eigen::VectorXd theta = _orthoG / _orthoG.norm();
      _previousOrthoG.noalias() = _orthoG;
      _previousOrthoFN.noalias() = _orthoFN;
      _oldParametersR1.noalias() = _parametersR1;
      /* Estimate curvature at phi = 0 from gradients of dimer */
      const double c0 = (_gradientsR1 - _gradients).dot(_dimerAxis) / radius;
      /* Estimate derivative of curvature relative to phi at phi = 0 from gradients of dimer */
      const double dCdphi = 2.0 * (_gradientsR1 - _gradients).dot(theta) / radius;
      /* check whether gradient at phi = 0 is already low enough to stop rotation */
      /* reason is that C at phi = 0 might already be lowest point, but because of too large phi1, this is not
       * recognized, therefore additional possibilty to stop rotation */
      if (fabs(dCdphi) < rotationGradientThresholdFirstCycle) {
        reachedTolerance = true;
        _sumRotations += rotationCycles;
        /* if already at first step and therefore no rotation is performed, the interval of rotations is increased */
        if (rotationCycles == 0)
          intervalOfRotations++;
        break;
      }
      /* estimate good phi for numerical minimum search */
      phi1 = -0.5 * std::atan(dCdphi / (2.0 * std::fabs(c0)));
      /* if already very low, rotation is converged (hardly the case) */
      if (fabs(phi1) < phiTolerance) {
        reachedTolerance = true;
        _sumRotations += rotationCycles;
        break;
      }
      /* dimer is rotated by phi1 */
      _dimerAxis.noalias() = (orthogonalDimerAxis * std::cos(phi1) + theta * std::sin(phi1)) * radius;
      _dimerAxis.normalize();
      _parametersR1.noalias() = parameters + radius * _dimerAxis;
      /* calculate gradient and then curvature at phi1 for finding phi minimum */
      function(_parametersR1, value1, gradientsR1phi1);
      const double cPhi1 = (gradientsR1phi1 - _gradients).dot(_dimerAxis) / radius;
      /* determine Fourier Coefficients a0, a1 and b1 */
      const double b1 = 0.5 * dCdphi;
      const double a1 = (c0 - cPhi1 + b1 * std::sin(2 * phi1)) / (1 - std::cos(2 * phi1));
      const double a0 = 2 * (c0 - a1);
      phiMin = 0.5 * std::atan(b1 / a1);
      /* Estimate curvature at phi minimum */
      cPhiMin = a0 / 2.0 + a1 * std::cos(2.0 * phiMin) + b1 * std::sin(2.0 * phiMin);
      /* If estimate is higher than both other curvature --> we found maximum, because of periodicity, we can add
       * pi/2 */
      if (cPhiMin > c0 && cPhiMin > cPhi1)
        phiMin += M_PI / 2;
      /* if only curvature at phi1 is lower, we rather pick phi1 for the minimum */
      else if (cPhi1 < cPhiMin)
        phiMin = phi1;
      /* same with phi = 0, this corresponds to converged rotation */
      else if (c0 < cPhiMin) {
        phiMin = 0.0;
        cPhiMin = c0;
      }
      /* rotate back by phi1, because phiMin from original position */
      _dimerAxis.noalias() = orthogonalDimerAxis;
      _parametersR1.noalias() = parameters + radius * _dimerAxis;
      if (fabs(phiMin) < phiTolerance) {
        /* already converged, no further rotation necessary */
        reachedTolerance = true;
        _sumRotations += rotationCycles;
        break;
      }
      /* rotate by phiMin */
      _dimerAxis.noalias() = (orthogonalDimerAxis * std::cos(phiMin) + theta * std::sin(phiMin)) * radius;
      _dimerAxis.normalize();
      _parametersR1.noalias() = parameters + radius * _dimerAxis;
    } // end of for loop for rotation
    if (!reachedTolerance) {
      std::cout << "WARNING: Did not reach tolerance in rotation, last angle in radiant was: " << phiMin << std::endl;
      std::cout << "If this angle is very high, optimization might not converge or converge to wrong saddle point "
                << std::endl;
      _sumRotations += maxRotationCycles;
      /* if break, R1 has already been calculated, but without break, calculation necessary. */
      function(_parametersR1, value1, _gradientsR1);
      determineDirectionOfDimerToBeMaximized(value, parameters, value1, cPhiMin, function);
      /* Calculate curvature anyway in case the dimer was not rotated */
      cPhiMin = (_gradientsR1 - _gradients).dot(_dimerAxis) / radius;
      return cPhiMin;
    }
    else {
      determineDirectionOfDimerToBeMaximized(value, parameters, value1, cPhiMin, function);
      /* No change necessary, the curvature is already known */
      return cPhiMin;
    }
  };

  /**
   * @brief Rotate the dimer directly based on the difference of the forces orthogonal to the dimer axis
   *        requires 1 calculation per rotation, only finds closest local minimum
   *
   * @param parameters         The parameters to be optimized
   * @param value              The value of the optimized point
   * @param maxRotationCycles  The number of maximum allowed rotations
   * @tparam UpdateFunction    A lambda function with a void return value, and the arguments:\n
   *                           1. const Eigen::VectorXd& parameters\n
   *                           2. double& value\n
   *                           3. Eigen::VectorXd& gradients
   *
   * @return double            Returns the curvature estimate of the rotated dimer
   */
  template<class UpdateFunction>
  double rotationWithGradient(const Eigen::VectorXd& parameters, const double& value, const int& maxRotationCycles,
                              UpdateFunction&& function) {
    Eigen::VectorXd preRotationDimerAxis = _dimerAxis;
    rotationGradientThreshold = rotationGradientThresholdOtherCycles;
    unsigned int m = 0; // for lbfgs
    double value1 = 0.0;
    double curvature = 0.0;
    /* option available to decrease the threshold for the rotation gradient after certain number of performed full
     * rotation cylces */
    if (decreaseRotationGradientThreshold && _numberOfPerformedRotationCycles > cycleOfRotationGradientDecrease)
      rotationGradientThreshold = loweredRotationGradientThreshold;
    /* starting individual rotations */
    for (int rotationCycles = 0; rotationCycles < maxRotationCycles; ++rotationCycles) {
      _orthoFN.noalias() =
          2.0 * ((_gradientsR1 - _gradients).dot(_dimerAxis)) * _dimerAxis - 2.0 * (_gradientsR1 - _gradients);
      _orthoG.noalias() = _orthoFN;
      /* Conjugate Gradient option for rotation */
      if (rotationCG)
        _orthoG = conjugateGradientRotation(rotationCycles);
      /* L-BFGS for rotation */
      else if (rotationLBFGS)
        _orthoG = lbfgsRotation(rotationCycles, m);
      /* Save old values */
      _previousOrthoG.noalias() = _orthoG;
      _previousOrthoFN.noalias() = _orthoFN;
      _oldParametersR1.noalias() = _parametersR1;
      /* Move R1, determine new dimer axis and R1 so the dimer length stays identical */
      _parametersR1 += _orthoG;
      _dimerAxis.noalias() = _parametersR1 - parameters;
      _dimerAxis.normalize();
      _parametersR1.noalias() = parameters + _dimerAxis * radius;
      function(_parametersR1, value1, _gradientsR1);
      curvature = (_gradientsR1 - _gradients).dot(_dimerAxis) / radius;
      if (_orthoG.norm() < rotationGradientThreshold) {
        _sumRotations += rotationCycles + 1;
        if (_dimerAxis.dot(preRotationDimerAxis) > 0.99)
          intervalOfRotations++;
        if (_dimerAxis.dot(preRotationDimerAxis) > 0.999)
          rotationGradientThreshold = rotationGradientThresholdOtherCycles * 10;
        else
          rotationGradientThreshold = rotationGradientThresholdOtherCycles;
        determineDirectionOfDimerToBeMaximized(value, parameters, value1, curvature, function);
        return curvature;
      }
    } // end of for loop for rotation
    std::cout << "WARNING: Did not reach tolerance in rotation, norm of last rotation gradient was: " << _orthoG.norm()
              << std::endl;
    std::cout << "If this is very high, the optimization might not converge or converge to wrong saddle point " << std::endl;
    intervalOfRotations = 1;
    _sumRotations += maxRotationCycles;
    determineDirectionOfDimerToBeMaximized(value, parameters, value1, curvature, function);
    return curvature;
  };

  /**
   * @brief Scale the step size based on projection change of the modified force onto the stepvector
   *
   *
   * @param parameters         The parameters to be optimized
   * @param value              The value of the optimized point
   * @tparam UpdateFunction    A lambda function with a void return value, and the arguments:\n
   *                           1. const Eigen::VectorXd& parameters\n
   *                           2. double& value\n
   *                           3. Eigen::VectorXd& gradients
   * @param previousMaxScaling The boolean whether maxscaling was used in the last scaling
   *
   * @return double            Returns the curvature estimate of the rotated dimer
   */
  template<class UpdateFunction>
  double scaleWithProjectionTrustRadius(Eigen::VectorXd& parameters, double& value, UpdateFunction&& function,
                                        bool& previousMaxScaling) {
    double currentStepLength = defaultTranslationStep;
    _modForce.noalias() = 2.0 * (_gradients.dot(_dimerAxis) * _dimerAxis) - _gradients;
    double projection1 = 1.0; // used for trust radius based on projection, always one, if stepvector not modified
    double projection0 = _modForce.dot(_G) / (_modForce.norm() * _G.norm());
    if (projection1 - projection0 < 1e-8)
      return defaultTranslationStep * 10;
    if (multiScale) {
      /* use with multiple scaleFactors */
      double maxScaling = 4; // empirical value
      double extrapolatedStepsizeFactor = -projection1 / (projection0 - projection1);
      if (projection0 > 0.0) {
        if (extrapolatedStepsizeFactor < maxScaling) {
          if (previousMaxScaling) {
            currentStepLength = extrapolatedStepsizeFactor * defaultTranslationStep;
            previousMaxScaling = false;
          }
          else
            currentStepLength *= extrapolatedStepsizeFactor;
        }
        else {
          currentStepLength *= maxScaling;
          previousMaxScaling = true;
        }
      }
      else if (projection0 < -0.5) {
        parameters += extrapolatedStepsizeFactor * _steps - _steps;
        function(parameters, value, _gradients);
      }
      else if (projection0 < 0.0)
        currentStepLength = extrapolatedStepsizeFactor * defaultTranslationStep;
    }
    else {
      /* always scale default stepsize */
      double extrapolatedStepsizeFactor = -projection1 / (projection0 - projection1);
      double maxScaling = 10; // empirical value, bigger here because no stacking of multiple factors
      if (projection0 > 0.0) {
        if (extrapolatedStepsizeFactor < maxScaling)
          currentStepLength = extrapolatedStepsizeFactor * defaultTranslationStep;
        else
          currentStepLength = maxScaling * defaultTranslationStep;
      }
      else if (projection0 < -0.5) {
        parameters += extrapolatedStepsizeFactor * _steps - _steps;
        function(parameters, value, _gradients);
      }
      else if (projection0 < 0.0)
        currentStepLength = extrapolatedStepsizeFactor * defaultTranslationStep;
    }
    return currentStepLength;
  };

  //// @brief The axis of the dimer
  Eigen::VectorXd _dimerAxis;
  //// @brief The gradients at the point to be optimized
  Eigen::VectorXd _gradients;
  //// @brief The parameters at R1
  Eigen::VectorXd _parametersR1;
  Eigen::VectorXd _oldParametersR1;
  //// @brief The gradients at R1
  Eigen::VectorXd _gradientsR1;
  //// @brief The projection of the force onto the dimer axis
  Eigen::VectorXd _fPara;
  Eigen::VectorXd _fParaOld;
  //// @brief L-BFGS matrix of parameter difference
  Eigen::MatrixXd _dx;
  //// @brief L-BFGS matrix of gradient difference
  Eigen::MatrixXd _dg;
  //// @brief The Force modified for TS search
  Eigen::VectorXd _modForce;
  //// @brief The stepvector of the translation
  Eigen::VectorXd _G;
  //// @brief The stepvector of rotation
  Eigen::VectorXd _orthoG;
  Eigen::VectorXd _previousOrthoG;
  //// @brief The difference of forces at the dimer points orthogonal to the dimer axis
  Eigen::VectorXd _orthoFN;
  Eigen::VectorXd _previousOrthoFN;
  //// @brief The scaled stepvector
  Eigen::VectorXd _steps;
  //// @brief The number of performed rotation cycles.
  unsigned int _numberOfPerformedRotationCycles = 0;
  //// @brief The number of performed single rotations.
  unsigned int _sumRotations = 0;
};
} // namespace Utils
} // namespace Scine

#endif // UTILS_DIMER_H_
