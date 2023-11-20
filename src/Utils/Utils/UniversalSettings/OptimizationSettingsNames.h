/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_OPTIMIZATIONSETTINGSNAMES_H
#define UTILS_OPTIMIZATIONSETTINGSNAMES_H

namespace Scine {
namespace Utils {
namespace SettingsNames {
namespace Optimizations {

//! SettingsNames provides a consistent list of settings names throughout the whole module.

struct Convergence {
  static constexpr const char* stepMaxCoeff = "convergence_step_max_coefficient";
  static constexpr const char* stepRMS = "convergence_step_rms";
  static constexpr const char* gradMaxCoeff = "convergence_gradient_max_coefficient";
  static constexpr const char* gradRMS = "convergence_gradient_rms";
  static constexpr const char* deltaValue = "convergence_delta_value";
  static constexpr const char* maxIter = "convergence_max_iterations";
  static constexpr const char* requirement = "convergence_requirement";
};

struct Bfgs {
  static constexpr const char* minIter = "bfgs_min_iterations";
  static constexpr const char* useTrustRadius = "bfgs_use_trust_radius";
  static constexpr const char* trustRadius = "bfgs_trust_radius";
  static constexpr const char* useGdiis = "bfgs_use_gdiis";
  static constexpr const char* gdiisMaxStore = "bfgs_gdiis_max_store";
};

struct Dimer {
  static constexpr const char* skipFirstRotation = "dimer_skip_first_rotation";
  static constexpr const char* decreaseRotationGradientThreshold = "dimer_decrease_rotation_gradient_threshold";
  static constexpr const char* gradientInterpolation = "dimer_gradient_interpolation";
  static constexpr const char* rotationCG = "dimer_rotation_conjugate_gradient";
  static constexpr const char* rotationLBFGS = "dimer_rotation_lbfgs";
  static constexpr const char* onlyOneRotation = "dimer_only_one_rotation";
  static constexpr const char* translation = "dimer_translation";
  static constexpr const char* multiScale = "dimer_multi_scale";
  static constexpr const char* radius = "dimer_radius";
  static constexpr const char* phiTolerance = "dimer_phi_tolerance";
  static constexpr const char* rotationGradientThresholdFirstCycle = "dimer_rotation_gradient_first";
  static constexpr const char* rotationGradientThresholdOtherCycles = "dimer_rotation_gradient_other";
  static constexpr const char* loweredRotationGradientThreshold = "dimer_lowered_rotation_gradient";
  static constexpr const char* gradRMSDthreshold = "dimer_grad_rmsd_threshold";
  static constexpr const char* trustRadius = "dimer_trust_radius";
  static constexpr const char* defaultTranslationStep = "dimer_default_translation_step";
  static constexpr const char* maxRotationsFirstCycle = "dimer_max_rotations_first_cycle";
  static constexpr const char* maxRotationsOtherCycles = "dimer_max_rotations_other_cycle";
  static constexpr const char* intervalOfRotations = "dimer_interval_of_rotations";
  static constexpr const char* cycleOfRotationGradientDecrease = "dimer_cycle_of_rotation_gradient_decrease";
  static constexpr const char* lbfgsMemory = "dimer_lbfgs_memory";
  static constexpr const char* bfgsStart = "dimer_bfgs_start";
  static constexpr const char* minimizationCycle = "dimer_minimization_cycle";
};

struct Lbfgs {
  static constexpr const char* maxm = "lbfgs_maxm";
  static constexpr const char* linesearch = "lbfgs_linesearch";
  static constexpr const char* c1 = "lbfgs_c1";
  static constexpr const char* c2 = "lbfgs_c2";
  static constexpr const char* stepLength = "lbfgs_step_length";
  static constexpr const char* useTrustRadius = "lbfgs_use_trust_radius";
  static constexpr const char* trustRadius = "lbfgs_trust_radius";
  static constexpr const char* maxBacktracking = "lbfgs_max_backtracking";
};

struct SteepestDescent {
  static constexpr const char* factor = "sd_factor";
  static constexpr const char* useTrustRadius = "sd_use_trust_radius";
  static constexpr const char* trustRadius = "sd_trust_radius";
  static constexpr const char* dynamicMultiplier = "sd_dynamic_multiplier";
};

struct Bofill {
  static constexpr const char* trustRadius = "bofill_trust_radius";
  static constexpr const char* hessianUpdate = "bofill_hessian_update";
  static constexpr const char* mode = "bofill_follow_mode";
};

struct EigenvectorFollowing {
  static constexpr const char* trustRadius = "ev_trust_radius";
  static constexpr const char* mode = "ev_follow_mode";
};

struct NewtonRaphson {
  static constexpr const char* trustRadius = "nr_trust_radius";
  static constexpr const char* svdThreshold = "nr_svd_threshold";
};

struct GeometryOptimizer {
  static constexpr const char* fixedAtoms = "geoopt_constrained_atoms";
  static constexpr const char* coordinateSystem = "geoopt_coordinate_system";
};

struct CellOptimizer {
  static constexpr const char* optimizeAngles = "cellopt_optimize_angles";
  static constexpr const char* optimizeA = "cellopt_optimize_a";
  static constexpr const char* optimizeB = "cellopt_optimize_b";
  static constexpr const char* optimizeC = "cellopt_optimize_c";
  static constexpr const char* geooptMaxIterations = "cellopt_geoopt_max_convergence_iterations";
  static constexpr const char* celloptMaxIterations = "cellopt_cellopt_max_convergence_iterations";
};

struct Nt {
  static constexpr const char* rHSList = "nt_rhs_list";
  static constexpr const char* lHSList = "nt_lhs_list";
  static constexpr const char* attractive = "nt_attractive";
  static constexpr const char* totalForceNorm = "nt_total_force_norm";
  static constexpr const char* maxIter = "convergence_max_iterations";
  static constexpr const char* repulsiveStop = "convergence_repulsive_stop";
  static constexpr const char* attractiveStop = "convergence_attractive_stop";
  static constexpr const char* sdFactor = "sd_factor";
  static constexpr const char* useMicroCycles = "nt_use_micro_cycles";
  static constexpr const char* fixedNumberOfMicroCycles = "nt_fixed_number_of_micro_cycles";
  static constexpr const char* numberOfMicroCycles = "nt_number_of_micro_cycles";
  static constexpr const char* filterPasses = "nt_filter_passes";
  static constexpr const char* extractionCriterion = "nt_extraction_criterion";
  static constexpr const char* coordinateSystem = "nt_coordinate_system";
  static constexpr const char* fixedAtoms = "nt_constrained_atoms";
  static constexpr const char* movableSide = "nt_movable_side";
  static constexpr const char* extractHighest = "highest_maximum";
  static constexpr const char* extractFirst = "first_maximum";
};

struct Nt2 {
  static constexpr const char* assList = "nt_associations";
  static constexpr const char* dissList = "nt_dissociations";
  static constexpr const char* totalForceNorm = "nt_total_force_norm";
  static constexpr const char* maxIter = "convergence_max_iterations";
  static constexpr const char* attractiveStop = "convergence_attractive_stop";
  static constexpr const char* sdFactor = "sd_factor";
  static constexpr const char* useMicroCycles = "nt_use_micro_cycles";
  static constexpr const char* fixedNumberOfMicroCycles = "nt_fixed_number_of_micro_cycles";
  static constexpr const char* numberOfMicroCycles = "nt_number_of_micro_cycles";
  static constexpr const char* filterPasses = "nt_filter_passes";
  static constexpr const char* extractionCriterion = "nt_extraction_criterion";
  static constexpr const char* coordinateSystem = "nt_coordinate_system";
  static constexpr const char* fixedAtoms = "nt_constrained_atoms";
  static constexpr const char* extractLastBeforeTarget = "last_maximum_before_first_target";
  static constexpr const char* extractHighest = "highest_maximum";
  static constexpr const char* extractFirst = "first_maximum";
};

struct Irc {
  static constexpr const char* initialStepSize = "irc_initial_step_size";
  static constexpr const char* coordinateSystem = "irc_coordinate_system";
};

struct Afir {
  static constexpr const char* rHSList = "afir_rhs_list";
  static constexpr const char* lHSList = "afir_lhs_list";
  static constexpr const char* weakForces = "afir_weak_forces";
  static constexpr const char* attractive = "afir_attractive";
  static constexpr const char* energyAllowance = "afir_energy_allowance";
  static constexpr const char* phaseIn = "afir_phase_in";
  static constexpr const char* coordinateSystem = "afir_coordinate_system";
};

struct MachineLearning {
  static constexpr const char* restartOptimization = "restart_optimization";
  static constexpr const char* numRestarts = "num_restarts";
  static constexpr const char* maxIterations = "max_iterations";
  static constexpr const char* maxLinesearch = "max_linesearch";
  static constexpr const char* convergenceTolerance = "convergence_tolerance";
  static constexpr const char* ftol = "linesearch_tolerance";
};

} // namespace Optimizations
} // namespace SettingsNames
} // namespace Utils
} // namespace Scine

#endif // UTILS_OPTIMIZATIONSETTINGSNAMES_H
