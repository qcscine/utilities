/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_DIMER_H_
#define UTILS_DIMER_H_

#include "Utils/Optimizer/GradientBased/Gdiis.h"
#include "Utils/Optimizer/GradientBased/GradientBasedCheck.h"
#include "Utils/Optimizer/Optimizer.h"
#include <Core/Log.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>

namespace Scine {
namespace Utils {

/**
 * @brief An implementation of the Dimer optimization algorithm for saddlepoints.
 *
 * transition state search without hessian by finding minimum curvature by gradient calculations of a dimer on PES
 * dimer is constructed by given vector or otherwise a random vector. Then dimer is rotated iteratively until it is
 * aligned with the lowest curvature mode, then translation is performed
 * this is repeated until converged
 *
 * This implementation is a mixture of the descriptions of following publications:
 *          Kaestner J. and Sherwood P., Journal of Chemical Physics, 2008, DOI: 10.1063/1.2815812
 *          Shang C. and Liu ZP., Journal of Chemical Theory and Computation, 2010, DOI: 10.1021/ct9005147
 *
 * General procedure:
 * The Dimer consists of the point to be optimized and a point along the dimer axis with a provided distance
 * ('dimer_radius'). Because there are two possible points on the axis, this point is always the one with the higher
 * value.
 * At each x optimization cycle given by 'dimer_interval_of_rotations' a check is performed if a rotation of the dimer
 * should be performed. This check is adapted from Shang et al. The rotation at the very beginning is performed with
 * the algorithm described by Kaestner et al. and is named 'rotation_with_phi'. In all other rotations the algorithm by
 * Shang et al. with the modification of using an L-BFGS instead of Broyden algorithm. The method of Shang to directly
 * optimize the gradient along the rotation axis converges faster, but converges to the next local minimum. At least
 * on simple test surfaces the rotation algorithm by Kaestner et al., which minimizes the rotation angle via an L-BFGS
 * minimization of the rotationvector, found the global minimum. Because the dimer can be generated from a random vector
 * at the beginning, we therefore first use the algorithm by Kaestner et al. and then later the fast algorithm by Shang
 * et al., because the next local minimum should still correspond to the global minimum.
 * Then, the actual optimizations is performed via a translation. Thereby, all modes are minimized and only the mode
 * along the dimer axis is maximized. The curvature along the dimer can be approximated as we know the gradient at the
 * two points of the dimer. If this curvature is still positive, no minimization is performed and only this mode is
 * maximized. This restriction is lifted after a given number of cycles ('dimer_minimization_cycle') as we saw an
 * improvement in needed optimization cycles. The algorithm of Shang et al. which further changes the weighting of the
 * optimization along the dimer axis and all other modes, is not applied as we have not yet performed a systematic study
 * to retrieve reliable thresholds, which are needed for this.
 * For the translation, further modifications are implemented, which are specified by the 'dimer_translation' option.
 * The default is a BFGS algorithm based on the already modified gradient based on the dimer axis. It is only activated
 * after a certain number of cycles specified by 'dimer_bfgs_start'. The several options are the reason for the several
 * if/else clauses after the rotation as either only the steplength or the stepvector are modified by the selected
 * algorithm. Furthermore, a trust radius ('dimer_trust_radius') is applied for all but the GDIIS, which already has
 * several sanity checks in his functions. After the translation, an oscillation check is performed. If positive, the
 * steplength is decreased and the dimer is rotated.
 * All further options not mentioned here, cause rather minor changes and can probably only be optimized by systematic
 * studies.
 *
 */
class Dimer : public Optimizer {
 public:
  static constexpr const char* dimerSkipFirstRotation = "dimer_skip_first_rotation";
  static constexpr const char* dimerDecreaseRotationGradientThreshold = "dimer_decrease_rotation_gradient_threshold";
  static constexpr const char* dimerGradientInterpolation = "dimer_gradient_interpolation";
  static constexpr const char* dimerRotationCG = "dimer_rotation_conjugate_gradient";
  static constexpr const char* dimerRotationLBFGS = "dimer_rotation_lbfgs";
  static constexpr const char* dimerOnlyOneRotation = "dimer_only_one_rotation";
  static constexpr const char* dimerTranslation = "dimer_translation";
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
  static constexpr const char* dimerLbfgsMemory = "dimer_lbfgs_memory";
  static constexpr const char* dimerBfgsStart = "dimer_bfgs_start";
  static constexpr const char* dimerMinimizationCycle = "dimer_minimization_cycle";
  /// @brief Default constructor.
  Dimer() = default;
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
  int optimize(Eigen::VectorXd& parameters, UpdateFunction&& function, ConvergenceCheckClass& check, Core::Log& log) {
    /* number of parameters treated */
    unsigned int nParams = parameters.size();
    if (nParams == 0) {
      throw EmptyOptimizerParametersException();
    }
    setVectorsToZero(nParams);
    bool previousMaxScaling = false;
    double value = 0.0;
    double currentStepLength = defaultTranslationStep;
    sanityCheck();
    _cycle = _startCycle;
    /* First calculation */
    function(parameters, value, _gradients);
    _gradientCalls++;
    this->triggerObservers(_cycle, value, parameters);
    // Init parameters and value stored in the convergence checker
    check.setParametersAndValue(parameters, value);
    createDimerAxis(nParams, function, parameters, value, log);
    bool stop = false;
    while (!stop) {
      _cycle++;
      /* Generate Dimer */
      _parametersR1.noalias() = parameters + radius * _dimerAxis;
      /* Force parallel to the dimer axis */
      _fParaOld.noalias() = _fPara;
      _fPara.noalias() = -_gradients.dot(_dimerAxis) * _dimerAxis;
      /* check whether rotation is performed in this cycle */
      unsigned int maxRotationCycles = (_cycle == _startCycle + 1) ? maxRotationsFirstCycle : maxRotationsOtherCycles;
      if (determineIfPerformRotation(_cycle, parameters, function)) {
        _numberOfPerformedRotationCycles++;
        if (_cycle == _startCycle + 1) {
          rotationWithPhi(parameters, value, maxRotationCycles, function, log);
        }
        else {
          rotationWithGradient(parameters, value, maxRotationCycles, function, log);
        }
        /* Reset scaling after rotation */
        currentStepLength = defaultTranslationStep;
      }
      if (_curvature > 0 && (_cycle - _startCycle) < static_cast<int>(minimizationCycle)) {
        /* no negative eigenvalue estimate --> still in convex region, trying to solely maximize lowest mode */
        _modGradient.noalias() = -_gradients.dot(_dimerAxis) * _dimerAxis;
      }
      else {
        _modGradient.noalias() = -2.0 * (_gradients.dot(_dimerAxis) * _dimerAxis) + _gradients;
      }
      if (_useTranslationBFGS) {
        _stepvector.noalias() = bfgsTranslation(parameters);
      }
      else if (_useAmsGrad) {
        _stepvector.noalias() = amsGrad(_cycle);
      }
      else {
        _stepvector.noalias() = -_modGradient;
      }

      /* if modification to step size, gradient is checked and scaling is activated below certain threshold */
      /* if scaling has been activated, no need to calculate gradientRMSD, otherwise calculate */
      if (_useProjectionLineSearch && (_appliedStepsizeScaling || gradientBelowThreshold())) {
        _appliedStepsizeScaling = true; // bool turned true if entered by gradient calculations
        /* projection method does not change direction and length has been changed at the end of last cycle */
        _steps = currentStepLength * _stepvector;
      }
      /* scaling cannot be applied yet, apply default scaling of step */
      /* this also applies for other translation methods */
      else {
        _steps = defaultTranslationStep * _stepvector;
      }

      /* Check trustradius and take a step */
      double maxVal = _steps.array().abs().maxCoeff();
      if (maxVal > trustRadius) {
        /* Scale step vector to trust radius */
        _steps *= trustRadius / maxVal;
        /* reset BFGS / stepsize scaling */
        invH = 0.5 * Eigen::MatrixXd::Identity(nParams, nParams);
        _appliedStepsizeScaling = false;
        currentStepLength = defaultTranslationStep;
      }
      constrainedAdd(parameters, _steps);
      /* Calculate energy and gradient of new position */
      function(parameters, value, _gradients);
      _gradientCalls++;
      this->triggerObservers(_cycle, value, parameters);
      /* Check convergence */
      stop = check.checkMaxIterations(_cycle) || check.checkConvergence(parameters, value, constrainGradient(_gradients));
      if (!stop && this->isOscillating(value)) {
        _oldParameters.noalias() = parameters;
        oscillationCorrection(_steps, parameters);
        function(parameters, value, _gradients);
        _cycle++;
        _gradientCalls++;
        this->triggerObservers(_cycle, value, parameters);
        defaultTranslationStep *= 0.95;
        currentStepLength = defaultTranslationStep;
        rotationWithGradient(parameters, value, maxRotationCycles, function, log);
      }
      else {
        defaultTranslationStep = 1.0;
        currentStepLength = defaultTranslationStep;
        if (_useProjectionLineSearch && _appliedStepsizeScaling && !stop) {
          currentStepLength = scaleWithProjectionLineSearch(parameters, value, function, previousMaxScaling);
        }
      }
    } // end of optimization while loop
    log.output << Core::Log::nl << "    Number of individual rotations: " << _sumRotations << Core::Log::nl;
    log.output << Core::Log::nl << "    Number of gradient calls: " << _gradientCalls << Core::Log::nl;
    return _cycle;
  }

  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include the Dimers's options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
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
    UniversalSettings::OptionListDescriptor translation(
        "Which algorithm shall be used for the translation of the dimer.");
    translation.addOption("bfgs");
    translation.addOption("linesearch");
    translation.addOption("amsgrad");
    translation.setDefaultOption("bfgs");
    collection.push_back(Dimer::dimerTranslation, translation);
    UniversalSettings::BoolDescriptor multi_scale(
        "If the projection step size factor shall be applied to previous step size.");
    multi_scale.setDefaultValue(multiScale);
    collection.push_back(Dimer::dimerMultiScale, multi_scale);
    UniversalSettings::DoubleDescriptor dimer_radius("Distance between the two images of the dimer.");
    dimer_radius.setMinimum(0.0);
    dimer_radius.setDefaultValue(radius);
    collection.push_back(Dimer::dimerRadius, dimer_radius);
    UniversalSettings::DoubleDescriptor phi_tolerance("The convergence criterion for the rotation angle.");
    phi_tolerance.setMinimum(0.0);
    phi_tolerance.setDefaultValue(phiTolerance);
    collection.push_back(Dimer::dimerPhiTolerance, phi_tolerance);
    UniversalSettings::DoubleDescriptor rotation_gradient_threshold_first_cycle(
        "The convergence criterion for the gradient of the rotation in the first rotation cycle.");
    rotation_gradient_threshold_first_cycle.setMinimum(0.0);
    rotation_gradient_threshold_first_cycle.setDefaultValue(rotationGradientThresholdFirstCycle);
    collection.push_back(Dimer::dimerRotationGradientThresholdFirstCycle, rotation_gradient_threshold_first_cycle);
    UniversalSettings::DoubleDescriptor rotation_gradient_threshold_other_cycles(
        "The convergence criterion for the gradient of the rotation in all but the first rotation cycle.");
    rotation_gradient_threshold_other_cycles.setMinimum(0.0);
    rotation_gradient_threshold_other_cycles.setDefaultValue(rotationGradientThresholdOtherCycles);
    collection.push_back(Dimer::dimerRotationGradientThresholdOtherCycles, rotation_gradient_threshold_other_cycles);
    UniversalSettings::DoubleDescriptor lowered_rotation_gradient_threshold(
        "The convergence criterion for the gradient of the rotation if the criterion is lowered after some cycles.");
    lowered_rotation_gradient_threshold.setMinimum(0.0);
    lowered_rotation_gradient_threshold.setDefaultValue(loweredRotationGradientThreshold);
    collection.push_back(Dimer::dimerLoweredRotationGradientThreshold, lowered_rotation_gradient_threshold);
    UniversalSettings::DoubleDescriptor grad_rmsd_threshold(
        "RMSD threshold of the gradient for starting the stepsize scaling.");
    grad_rmsd_threshold.setMinimum(0.0);
    grad_rmsd_threshold.setDefaultValue(gradRMSDthreshold);
    collection.push_back(Dimer::dimerGradRMSDthreshold, grad_rmsd_threshold);
    UniversalSettings::DoubleDescriptor dimer_trust_radius("The maximum RMS of a taken step.");
    dimer_trust_radius.setMinimum(0.0);
    dimer_trust_radius.setDefaultValue(trustRadius);
    collection.push_back(Dimer::dimerTrustRadius, dimer_trust_radius);
    UniversalSettings::DoubleDescriptor dimer_default_translation_step(
        "The factor to multiple the stepsize vector in a steepest descent step.");
    dimer_default_translation_step.setMinimum(0.0);
    dimer_default_translation_step.setDefaultValue(defaultTranslationStep);
    collection.push_back(Dimer::dimerDefaultTranslationStep, dimer_default_translation_step);
    UniversalSettings::IntDescriptor max_rotations_first_cycle(
        "The maximum number of rotations in the first rotation cycle.");
    max_rotations_first_cycle.setMinimum(0);
    max_rotations_first_cycle.setDefaultValue(maxRotationsFirstCycle);
    collection.push_back(Dimer::dimerMaxRotationsFirstCycle, max_rotations_first_cycle);
    UniversalSettings::IntDescriptor max_rotations_other_cycle(
        "The maximum number of rotations in all but the first rotation cycle.");
    max_rotations_other_cycle.setMinimum(0);
    max_rotations_other_cycle.setDefaultValue(maxRotationsOtherCycles);
    collection.push_back(Dimer::dimerMaxRotationsOtherCycles, max_rotations_other_cycle);
    UniversalSettings::IntDescriptor interval_of_rotations(
        "The interval of performed rotation cycles in the total optimization steps.");
    interval_of_rotations.setMinimum(1);
    interval_of_rotations.setDefaultValue(intervalOfRotations);
    collection.push_back(Dimer::dimerIntervalOfRotations, interval_of_rotations);
    UniversalSettings::IntDescriptor cycle_of_rotation_gradient_decrease(
        "The number of rotation cycles after which the threshold for the gradient of rotation is decreased.");
    cycle_of_rotation_gradient_decrease.setMinimum(0);
    cycle_of_rotation_gradient_decrease.setDefaultValue(cycleOfRotationGradientDecrease);
    collection.push_back(Dimer::dimerCycleOfRotationGradientDecrease, cycle_of_rotation_gradient_decrease);
    UniversalSettings::IntDescriptor lbfgs_memory("The number of saved gradients during rotation in L-BFGS.");
    lbfgs_memory.setMinimum(1);
    lbfgs_memory.setDefaultValue(lbfgsMemory);
    collection.push_back(Dimer::dimerLbfgsMemory, lbfgs_memory);
    UniversalSettings::IntDescriptor bfgs_start("The cycle in which BFGS is used in translation.");
    bfgs_start.setMinimum(1);
    bfgs_start.setDefaultValue(bfgsStart);
    collection.push_back(Dimer::dimerBfgsStart, bfgs_start);
    UniversalSettings::IntDescriptor minimization_cycle("The cycle in which all other modes are always minimized");
    minimization_cycle.setMinimum(1);
    minimization_cycle.setDefaultValue(minimizationCycle);
    collection.push_back(Dimer::dimerMinimizationCycle, minimization_cycle);
  }
  /**
   * @brief Updates the Dimer's options with those values given in the Settings.
   * @param settings The settings to update the option of the steepest descent with.
   */
  void applySettings(const Settings& settings) final {
    skipFirstRotation = settings.getBool(Dimer::dimerSkipFirstRotation);
    decreaseRotationGradientThreshold = settings.getBool(Dimer::dimerDecreaseRotationGradientThreshold);
    gradientInterpolation = settings.getBool(Dimer::dimerGradientInterpolation);
    rotationCG = settings.getBool(Dimer::dimerRotationCG);
    rotationLBFGS = settings.getBool(Dimer::dimerRotationLBFGS);
    onlyOneRotation = settings.getBool(Dimer::dimerOnlyOneRotation);
    translationMethod = settings.getString(Dimer::dimerTranslation);
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
    lbfgsMemory = settings.getInt(Dimer::dimerLbfgsMemory);
    bfgsStart = settings.getInt(Dimer::dimerBfgsStart);
    minimizationCycle = settings.getInt(Dimer::dimerMinimizationCycle);
  }

  /**
   * @brief Prepares the Dimer optimizer for rerunning its optimize function.
   *
   * This function is used to prepare a Dimer optimizer instance for rerunning
   * its optimize function with possibly different settings.
   * It changes the cycle count the optimizer starts with when the optimize
   * function is called and removes the stored inverse Hessian and its
   * projection function. The optional guess vector used in the
   * first step is removed.
   * The value memory used for oscillating correction is cleared.
   *
   * @param cycleNumber The cycle number the optimizer starts with.
   */
  void prepareRestart(const int cycleNumber) final {
    // Set start cycle number
    _startCycle = cycleNumber;
    // Remove inverse Hessian and projection function
    (this->invH).resize(0, 0);
    this->projection = nullptr;
    // Remove optional starting guess vector
    this->guessVector = nullptr;
    // Clear value memory for oscillating correction
    _valueMemory.clear();
  }
  //// @brief If the provided guessVector shall be used without rotation in first step.
  bool skipFirstRotation = false;
  //// @brief If threshold for convergence in rotation shall be lowered after certain number of rotationcycles.
  bool decreaseRotationGradientThreshold = false;
  //// @brief If gradient of R1 shall be estimated by interpolation from last rotation data.
  bool gradientInterpolation = true;
  //// @brief If conjugate gradient shall be used for rotation.
  bool rotationCG = false;
  //// @brief If L-BFGS shall be used for rotation.
  bool rotationLBFGS = true;
  //// @brief If rotation shall only be performed in first step.
  bool onlyOneRotation = false;
  //// @brief If the projection step size factor shall be applied to previous step size.
  bool multiScale = true;
  //// @brief Which algorithm shall be used for translation. Possible values are BFGS, Linesearch, AMSGRAD
  std::string translationMethod = "BFGS";
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
  //// @brief The convergence criterion for the gradient of the rotation if the criterion is lowered after some cycles.
  double gradRMSDthreshold = 1e-3;
  /// @brief The maximum RMS of a taken step.
  double trustRadius = 0.2;
  /// @brief The factor to multiple the stepsize vector in a steepest descent step.
  double defaultTranslationStep = 1;
  //// @brief parameter for AMSGRAD
  double alpha1 = 0.01;
  //// @brief parameter for AMSGRAD
  double beta1 = 0.1;
  //// @brief parameter for AMSGRAD
  double beta2 = 0.01;
  //// @brief The maximum number of rotations in the first rotation cycle.
  unsigned int maxRotationsFirstCycle = 100;
  //// @brief The maximum number of rotations in all but the first rotation cycle.
  unsigned int maxRotationsOtherCycles = 100;
  //// @brief The interval of performed rotation cycles in the total optimization steps.
  unsigned int intervalOfRotations = 5;
  //// @brief The number of rotation cycles after which the threshold for the gradient of rotation is decreased.
  unsigned int cycleOfRotationGradientDecrease = 5;
  //// @brief The number of saved gradients during rotation in L-BFGS.
  unsigned int lbfgsMemory = 5;
  //// @brief The cycle in which BFGS is started for the translation.
  unsigned int bfgsStart = 16;
  //// @brief The cycle in which the minimization of all but the dimer direction is enforced
  unsigned int minimizationCycle = 5;
  //// @brief The (optional) dimer axis for the first step provided by hessian calculation or trajectory.
  std::shared_ptr<Eigen::VectorXd> guessVector;
  /// @brief A possible Hessian projection
  std::shared_ptr<std::function<void(Eigen::MatrixXd&)>> projection = nullptr;
  //// @brief inverse hessian matrix
  Eigen::MatrixXd invH;

 private:
  void setVectorsToZero(const int& nParams) {
    _learningVector = Eigen::VectorXd::Zero(nParams);
    _oldLearningVector = Eigen::VectorXd::Zero(nParams);
    _gradModification = Eigen::VectorXd::Zero(nParams);
    _oldGradients = Eigen::VectorXd::Zero(nParams);
    _oldParameters = Eigen::VectorXd::Zero(nParams);
    _gradients = Eigen::VectorXd::Zero(nParams);
    _gradientsR1 = Eigen::VectorXd::Zero(nParams);
    _fPara = Eigen::VectorXd::Zero(nParams);
    _fParaOld = Eigen::VectorXd::Zero(nParams);
    _modGradient = Eigen::VectorXd::Zero(nParams);
    _stepvector = Eigen::VectorXd::Zero(nParams);
    _orthoFN = Eigen::VectorXd::Zero(nParams);
    _orthoG = Eigen::VectorXd::Zero(nParams);
    _previousOrthoFN = Eigen::VectorXd::Zero(nParams);
    _previousOrthoG = Eigen::VectorXd::Zero(nParams);
    _steps = Eigen::VectorXd::Zero(nParams);
    _dx.resize(nParams, lbfgsMemory);
    _dg.resize(nParams, lbfgsMemory);
    _dx.setZero();
    _dg.setZero();
    /* Set inverse hessian to identity matrix if nothing has been provided */
    if (invH.size() == 0) {
      invH = 0.5 * Eigen::MatrixXd::Identity(nParams, nParams);
    }
  }

  void sanityCheck() {
    // change input string to lowercase
    std::for_each(translationMethod.begin(), translationMethod.end(), [](char& c) { c = ::tolower(c); });
    if (translationMethod == "bfgs") {
      _useTranslationBFGS = true;
    }
    else if (translationMethod == "linesearch") {
      _useProjectionLineSearch = true;
    }
    else if (translationMethod == "amsgrad") {
      _useAmsGrad = true;
    }
    else {
      throw std::runtime_error(
          "Option " + translationMethod +
          " is not known for the translation of the dimer. Either use 'bfgs', 'linesearch', or 'amsgrad'.");
    }

    if (rotationLBFGS && rotationCG) {
      /* LBFGS is default, so CG is wished to be performed */
      rotationLBFGS = false;
    }
  }

  /**
   * @brief Create first dimer
   *
   * If guessvector is given, it is used for the dimer axis
   * if no guess, a random vector is used
   *
   * @param nParams            The number of parameters to be optimized
   * @tparam UpdateFunction    A lambda function with a void return value, and the arguments:\n
   *                           1. const Eigen::VectorXd& parameters\n
   *                           2. double& value\n
   *                           3. Eigen::VectorXd& gradients
   * @param parameters         The parameters to be optimized
   * @param value              The value of the optimized point
   *
   * @return void
   */
  template<class UpdateFunction>
  void createDimerAxis(const int& nParams, UpdateFunction&& function, const Eigen::VectorXd& parameters,
                       const double value, Core::Log& log) {
    bool randomNecessary = false;
    if (guessVector) {
      if ((*guessVector).size() == nParams) {
        _dimerAxis.noalias() = *guessVector;
      }
      else {
        log.warning << "Warning, given guess vector has different size than parameters. "
                       "Continuing transition state search, but omitting the guess."
                    << Core::Log::nl;
        guessVector = nullptr;
        randomNecessary = true;
      }
    }
    else {
      randomNecessary = true;
    }
    if (randomNecessary) {
      srand(42);
      _dimerAxis.noalias() = Eigen::VectorXd::Random(nParams);
    }
    _dimerAxis.normalize();
    /* The given vector may not be aligned in the right direction, but the opposite one */
    /* Therefore, construct R1, check value, then invert axis if value decreases */
    _parametersR1.noalias() = parameters + radius * _dimerAxis;
    double valueR1 = 0.0;
    function(_parametersR1, valueR1, _gradientsR1);
    _gradientCalls++;
    determineDirectionOfDimerToBeMaximized(value, parameters, valueR1, function);
  }

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
  bool determineIfPerformRotation(const int& cycle, const Eigen::VectorXd& /* parameters */, UpdateFunction&& function) {
    double value1;
    if (cycle == _startCycle + 1) {
      if (skipFirstRotation) {
        return false;
      }

      _rotationGradientThreshold = rotationGradientThresholdFirstCycle;
      function(_parametersR1, value1, _gradientsR1);
      _gradientCalls++;
      /* determine ortho force at current position to rotate */
      _orthoFN = 2.0 * ((_gradientsR1 - _gradients).dot(_dimerAxis)) * _dimerAxis - 2.0 * (_gradientsR1 - _gradients);
      _curvature = (_gradientsR1 - _gradients).dot(_dimerAxis) / radius;
      return true;
    }
    /* if only first rotation wants to be performed */
    if (onlyOneRotation) {
      return false;
    }

    if ((cycle - _startCycle) % intervalOfRotations == 0) {
      _rotationGradientThreshold = rotationGradientThresholdOtherCycles;
      function(_parametersR1, value1, _gradientsR1);
      _gradientCalls++;
      _curvature = (_gradientsR1 - _gradients).dot(_dimerAxis) / radius;
      /* determine ortho force at current position to rotate */
      _orthoFN = 2.0 * ((_gradientsR1 - _gradients).dot(_dimerAxis)) * _dimerAxis - 2.0 * (_gradientsR1 - _gradients);
      if ((_curvature > 0 && _fPara.norm() < _fParaOld.norm()) || (_curvature < 0 && _fPara.norm() > _fParaOld.norm())) {
        return true;
      }

      if (_orthoFN.norm() > 1e-3) {
        return true;
      }
    }
    return false;
  }

  /**
   * @brief Determine step vector of rotation with the Conjugate Gradient (Polak Ribiere) method
   *
   * @param rotationCycles     The number of individual rotation in the current rotation - 1
   *
   * @return Eigen::VectorXd   Returns step vector of rotation
   */
  Eigen::VectorXd conjugateGradientRotation(const unsigned rotationCycles) {
    _orthoG = _orthoFN;
    if (rotationCycles != 0) {
      /* Polak Ribiere */
      double beta = _orthoFN.dot(_orthoFN - _previousOrthoFN) / (_previousOrthoFN.dot(_previousOrthoFN));
      _orthoG = _orthoFN + beta * _previousOrthoG;
    }
    return _orthoG;
  }

  /**
   * @brief Determine the step vector of rotation based on the L-BFGS methods
   *
   * @param rotationCycles     The number of individual rotation in the current rotation - 1
   * @param m                  The counter for the L-BFGS backtracking
   *
   * @return Eigen::VectorXd   Returns step vector of rotation
   */
  Eigen::VectorXd lbfgsRotation(const unsigned rotationCycles, unsigned int& m) {
    _orthoG = _orthoFN;
    unsigned maxm = lbfgsMemory;
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
    for (int i = static_cast<int>(m) - 1; i > -1; --i) {
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
  }

  /**
   * @brief Flips dimeraxis if energy decreases along axis
   *
   * Rotate dimer axis if necessary, so dimer axis is afterwards correctly aligned with the direction in which the mode
   * can be optimized
   *
   * @param value              The value of the optimized point
   * @param parameters         The parameters to be optimized
   * @param value1             The value at R1
   * @tparam UpdateFunction    A lambda function with a void return value, and the arguments:\n
   *                           1. const Eigen::VectorXd& parameters\n
   *                           2. double& value\n
   *                           3. Eigen::VectorXd& gradients
   *
   * @return void
   */
  template<class UpdateFunction>
  void determineDirectionOfDimerToBeMaximized(const double& value, const Eigen::VectorXd& parameters, double& value1,
                                              UpdateFunction&& function) {
    /* if value1 bigger than value, do nothing, otherwise: */
    if (value > value1) {
      double curvatureOld = (_gradientsR1 - _gradients).dot(_dimerAxis) / radius;
      /* flip axis, calculate new point */
      _dimerAxis *= -1;
      Eigen::VectorXd parametersR1New = parameters + _dimerAxis * radius;
      Eigen::VectorXd gradientsR1New = _gradientsR1;
      double value1New = 0.0;
      function(parametersR1New, value1New, gradientsR1New);
      _gradientCalls++;
      /* if new point also lower in value, check curvature of both possible directions */
      if (value > value1New) {
        double curvatureNew = (gradientsR1New - _gradients).dot(_dimerAxis) / radius;
        /* go with lower curvature */
        if (curvatureNew < curvatureOld) {
          /* flipped axis has lower curvature, save values into class and exit */
          value1 = value1New;
          _parametersR1.noalias() = parametersR1New;
          _gradientsR1.noalias() = gradientsR1New;
          _curvature = curvatureNew;
          return;
        }

        /* previous direction had lower curvature -> flip back and exit */
        _dimerAxis *= -1;
        return;
      }

      /* new flipped value was higher in energy, while the old one wasn't -> save values and exit */
      value1 = value1New;
      _parametersR1.noalias() = parametersR1New;
      _gradientsR1.noalias() = gradientsR1New;
      _curvature = (_gradientsR1 - _gradients).dot(_dimerAxis) / radius;
    }
  }

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
   * @return void
   */
  template<class UpdateFunction>
  void rotationWithPhi(const Eigen::VectorXd& parameters, const double value, const unsigned int maxRotationCycles,
                       UpdateFunction&& function, Core::Log& log) {
    bool reachedTolerance = false;
    unsigned int m = 0; // for lbfgs
    double cPhiMin = 0.0;
    double phiMin = 0.0;
    double phi1 = 0.0;
    double value1 = 0.0;
    Eigen::VectorXd gradientsR1phi1 = Eigen::VectorXd::Zero(parameters.size());
    /* starting individual rotations */
    for (unsigned rotationCycles = 0; rotationCycles < maxRotationCycles; ++rotationCycles) {
      /* calculate E and F at R1 */
      /* either estimate gradient by interpolation or calculate it */
      if (gradientInterpolation && rotationCycles != 0) {
        _gradientsR1 = std::sin(phi1 - phiMin) / std::sin(phi1) * _gradientsR1 +
                       std::sin(phiMin) / std::sin(phi1) * gradientsR1phi1 +
                       (1 - std::cos(phiMin) - std::sin(phiMin) * std::tan(phi1 / 2)) * _gradients;
      }
      else {
        function(_parametersR1, value1, _gradientsR1);
        _gradientCalls++;
      }
      /* determine ortho force at current position to determin rotation plane */
      Eigen::VectorXd orthogonalDimerAxis = _dimerAxis; // used for transformation of the dimer axis
      _orthoFN.noalias() =
          2.0 * ((_gradientsR1 - _gradients).dot(_dimerAxis)) * _dimerAxis - 2.0 * (_gradientsR1 - _gradients);
      _orthoG.noalias() = _orthoFN;
      /* Conjugate Gradient option for rotation */
      if (rotationCG) {
        _orthoG = conjugateGradientRotation(rotationCycles);
        /* L-BFGS for rotation */
      }
      else if (rotationLBFGS) {
        _orthoG = lbfgsRotation(rotationCycles, m);
      }
      _orthoG -= _orthoG.dot(_dimerAxis) * _dimerAxis; // ensures that modified direction is still orthogonal to dimer
                                                       // axis
      /* unit vector which spans rotation plane together with dimer axis */
      Eigen::VectorXd theta = _orthoG / _orthoG.norm();
      if (theta != theta) {
        throw std::runtime_error("Gradient of zero length in Dimer rotation.");
      }
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
        if (rotationCycles == 0) {
          intervalOfRotations++;
        }
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
      _gradientCalls++;
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
      if (cPhiMin > c0 && cPhiMin > cPhi1) {
        phiMin += M_PI / 2;
        _curvature = a0 / 2.0 + a1 * std::cos(2.0 * phiMin) + b1 * std::sin(2.0 * phiMin);
      }
      /* if only curvature at phi1 is lower, we rather pick phi1 for the minimum */
      else if (cPhi1 < cPhiMin) {
        phiMin = phi1;
        _curvature = cPhi1;
      }
      /* same with phi = 0, this corresponds to converged rotation */
      else if (c0 < cPhiMin) {
        phiMin = 0.0;
        _curvature = c0;
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
      log.warning << "Warning, did not reach tolerance in rotation, last angle in radiant was: " << phiMin << Core::Log::nl
                  << "If this angle is very high, optimization might not converge or converge to wrong saddle point "
                  << Core::Log::nl;
      _sumRotations += maxRotationCycles;
      /* if break, R1 has already been calculated, but without break, calculation necessary. */
      function(_parametersR1, value1, _gradientsR1);
      _gradientCalls++;
      determineDirectionOfDimerToBeMaximized(value, parameters, value1, function);
      /* Calculate curvature anyway in case the dimer was not rotated */
      _curvature = (_gradientsR1 - _gradients).dot(_dimerAxis) / radius;
    }
    else {
      determineDirectionOfDimerToBeMaximized(value, parameters, value1, function);
      /* No calculation necessary, the curvature is already known */
    }
  }

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
   * @return void
   */
  template<class UpdateFunction>
  void rotationWithGradient(const Eigen::VectorXd& parameters, const double value, const unsigned int maxRotationCycles,
                            UpdateFunction&& function, Core::Log& log) {
    Eigen::VectorXd preRotationDimerAxis = _dimerAxis;
    _rotationGradientThreshold = rotationGradientThresholdOtherCycles;
    unsigned int m = 0; // for lbfgs
    double value1 = 0.0;
    /* option available to decrease the threshold for the rotation gradient after certain number of performed full
     * rotation cycles */
    if (decreaseRotationGradientThreshold && _numberOfPerformedRotationCycles > cycleOfRotationGradientDecrease) {
      _rotationGradientThreshold = loweredRotationGradientThreshold;
    }
    /* starting individual rotations */
    for (unsigned rotationCycles = 0; rotationCycles < maxRotationCycles; ++rotationCycles) {
      _orthoFN.noalias() =
          2.0 * ((_gradientsR1 - _gradients).dot(_dimerAxis)) * _dimerAxis - 2.0 * (_gradientsR1 - _gradients);
      _orthoG.noalias() = _orthoFN;
      /* Conjugate Gradient option for rotation */
      if (rotationCG) {
        _orthoG = conjugateGradientRotation(rotationCycles);
        /* L-BFGS for rotation */
      }
      else if (rotationLBFGS) {
        _orthoG = lbfgsRotation(rotationCycles, m);
      }
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
      _gradientCalls++;
      _curvature = (_gradientsR1 - _gradients).dot(_dimerAxis) / radius;
      if (_orthoG.norm() < _rotationGradientThreshold) {
        _sumRotations += rotationCycles + 1;
        if (_dimerAxis.dot(preRotationDimerAxis) > 0.99) {
          intervalOfRotations++;
        }
        if (_dimerAxis.dot(preRotationDimerAxis) > 0.999) {
          _rotationGradientThreshold = rotationGradientThresholdOtherCycles * 10;
        }
        else {
          _rotationGradientThreshold = rotationGradientThresholdOtherCycles;
        }
        determineDirectionOfDimerToBeMaximized(value, parameters, value1, function);
        return;
      }
    } // end of for loop for rotation
    log.warning << "Warning, did not reach tolerance in rotation, norm of last rotation gradient was: " << _orthoG.norm()
                << Core::Log::nl << "If this is very high, the optimization might not converge or converge to wrong saddle point "
                << Core::Log::nl;
    intervalOfRotations = 1;
    _sumRotations += maxRotationCycles;
    determineDirectionOfDimerToBeMaximized(value, parameters, value1, function);
  }

  bool gradientBelowThreshold() {
    return (std::sqrt(_gradients.squaredNorm() / _gradients.size()) < gradRMSDthreshold);
  }

  /**
   * @brief Determine stepvector of translation by BFGS algorithm
   *
   * @param parameters         The parameters to be optimized
   *
   * @return Eigen::VectorXd   Returns the stepvector
   */
  Eigen::VectorXd bfgsTranslation(const Eigen::VectorXd& parameters) {
    /* if first cycle save values and return unmodified vector */
    if (_cycle - _startCycle == 1) {
      _oldParameters.noalias() = parameters;
      _oldGradients.noalias() = _modGradient;
      return -_modGradient;
    }
    /* BFGS inverse Hessian update */
    Eigen::VectorXd dx = parameters - _oldParameters;
    Eigen::VectorXd dg = _modGradient - _oldGradients;
    double dxTdg = dx.dot(dg);
    const double dgTinvHdg = dg.transpose() * invH * dg;
    if (invH.isApprox(0.5 * Eigen::MatrixXd::Identity(parameters.size(), parameters.size()))) {
      /* Set initial inverse Hessian to dxTdg / dgTdg * I */
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
    if (std::fabs(delta - 1.0) < 1e-16) {
      dx.noalias() = delta * dx + (1.0 - delta) * invH * dg;
      dxTdg = dx.dot(dg);
    }

    /* Sanity check for dxTdg */
    if (fabs(dxTdg) < 1e-9) {
      dxTdg = (dxTdg < 0.0) ? -1e-9 : 1e-9;
    }
    const double alpha = (dxTdg + dg.transpose() * invH * dg) / (dxTdg * dxTdg);
    const double beta = 1.0 / dxTdg;
    invH += alpha * (dx * dx.transpose()) - beta * (invH * dg * dx.transpose() + dx * dg.transpose() * invH);

    /* save old values */
    _oldParameters.noalias() = parameters;
    _oldGradients.noalias() = _modGradient;
    if (projection) {
      (*projection)(invH);
    }
    if (!_appliedStepsizeScaling && (_cycle - _startCycle >= static_cast<int>(bfgsStart) || gradientBelowThreshold())) {
      _appliedStepsizeScaling = true;
    }
    if (_appliedStepsizeScaling) {
      return -invH * _modGradient;
    }
    return -_modGradient;
  }

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
   * @return double            Returns the step length
   */
  template<class UpdateFunction>
  double scaleWithProjectionLineSearch(Eigen::VectorXd& parameters, double& value, UpdateFunction&& function,
                                       bool& previousMaxScaling) {
    double currentStepLength = defaultTranslationStep;
    _modGradient.noalias() = -2.0 * (_gradients.dot(_dimerAxis) * _dimerAxis) + _gradients;
    double projection0 = -_modGradient.dot(_stepvector) / (_modGradient.norm() * _stepvector.norm());
    double projection1 = 1.0;
    if (projection1 - projection0 < 1e-8) {
      return defaultTranslationStep * 10;
    }
    if (multiScale) {
      /* use with multiple scaleFactors */
      double maxScaling = 2; // empirical value taken from Kästner paper
      double extrapolatedStepsizeFactor = -projection1 / (projection0 - projection1);
      if (projection0 > 0.0) {
        if (extrapolatedStepsizeFactor < maxScaling) {
          if (previousMaxScaling) {
            currentStepLength = extrapolatedStepsizeFactor * defaultTranslationStep;
            previousMaxScaling = false;
          }
          else {
            currentStepLength *= extrapolatedStepsizeFactor;
          }
        }
        else {
          currentStepLength *= maxScaling;
          previousMaxScaling = true;
        }
      }
      else if (projection0 < -0.5) {
        constrainedAdd(parameters, extrapolatedStepsizeFactor * _steps - _steps);
        function(parameters, value, _gradients);
        _gradientCalls++;
      }
      else if (projection0 < 0.0) {
        currentStepLength = extrapolatedStepsizeFactor * defaultTranslationStep;
      }
    }
    else {
      /* always scale default stepsize */
      double extrapolatedStepsizeFactor = -projection1 / (projection0 - projection1);
      double maxScaling = 10; // empirical value, bigger here because no stacking of multiple factors
      if (projection0 > 0.0) {
        if (extrapolatedStepsizeFactor < maxScaling) {
          currentStepLength = extrapolatedStepsizeFactor * defaultTranslationStep;
        }
        else {
          currentStepLength = maxScaling * defaultTranslationStep;
        }
      }
      else if (projection0 < -0.5) {
        constrainedAdd(parameters, extrapolatedStepsizeFactor * _steps - _steps);
        function(parameters, value, _gradients);
        _gradientCalls++;
      }
      else if (projection0 < 0.0) {
        currentStepLength = extrapolatedStepsizeFactor * defaultTranslationStep;
      }
    }
    return currentStepLength;
  }

  Eigen::VectorXd amsGrad(const int cycle) {
    _gradModification *= (1 - beta1);
    _gradModification += beta1 * _modGradient;
    _learningVector = (1 - beta2) * _oldLearningVector;
    for (int i = 0; i < _learningVector.size(); ++i) {
      _learningVector[i] += beta2 * _stepvector[i] * _stepvector[i];
    }
    _learningVector.cwiseMax(_oldLearningVector);
    _oldLearningVector = _learningVector;
    /* now if no scaling save values and return unmodified vector */
    if (cycle == _startCycle + 1) {
      return -_modGradient;
    }

    Eigen::MatrixXd learningMatrix = _learningVector.asDiagonal();
    learningMatrix = learningMatrix.cwiseSqrt();
    return -alpha1 * learningMatrix.inverse() * _gradModification;
  }
  /// @brief The counter variable for optimisation cycles
  int _cycle = 0;
  //// @brief The estimated curvature of the PES along the dimer
  double _curvature = 0.0;
  //// @brief The employed threshold depending on conditions.
  double _rotationGradientThreshold = 0.0;
  //// @brief The axis of the dimer
  Eigen::VectorXd _dimerAxis;
  //// @brief The gradients at the point to be optimized
  Eigen::VectorXd _gradients;
  Eigen::VectorXd _oldGradients;
  //// @brief The parameters at R1
  Eigen::VectorXd _parametersR1;
  Eigen::VectorXd _oldParametersR1;
  Eigen::VectorXd _oldParameters;
  //// @brief The gradients at R1
  Eigen::VectorXd _gradientsR1;
  //// @brief The projection of the force onto the dimer axis
  Eigen::VectorXd _fPara;
  Eigen::VectorXd _fParaOld;
  //// @brief L-BFGS matrix of parameter difference
  Eigen::MatrixXd _dx;
  //// @brief L-BFGS matrix of gradient difference
  Eigen::MatrixXd _dg;
  //// @brief The gradient modified for TS search
  Eigen::VectorXd _modGradient;
  //// @brief The stepvector of the translation
  Eigen::VectorXd _stepvector;
  //// @brief The stepvector of rotation
  Eigen::VectorXd _orthoG;
  Eigen::VectorXd _previousOrthoG;
  //// @brief The difference of forces at the dimer points orthogonal to the dimer axis
  Eigen::VectorXd _orthoFN;
  Eigen::VectorXd _previousOrthoFN;
  //// @brief The scaled stepvector
  Eigen::VectorXd _steps;
  Eigen::VectorXd _gradModification;
  Eigen::VectorXd _learningVector;
  Eigen::VectorXd _oldLearningVector;
  //// @brief The number of performed rotation cycles.
  unsigned int _numberOfPerformedRotationCycles = 0;
  //// @brief The number of performed single rotations.
  unsigned int _sumRotations = 0;
  unsigned int _gradientCalls = 0;
  bool _appliedStepsizeScaling = false;
  //// @brief If the step size shall be estimated by projection of the gradient at current step onto the last dimer axis.
  bool _useProjectionLineSearch = false;
  //// @brief If the step vector shall be modified by the AMSGRAD algorithm
  bool _useAmsGrad = false;
  //// @brief If the step vector shall be modified by BFGS algorithm
  bool _useTranslationBFGS = false;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_DIMER_H_
