/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Pybind.h>
#include <Utils/UniversalSettings/OptimizationSettingsNames.h>

using namespace Scine::Utils;

void init_opt_settings_names(pybind11::module& m) {
  pybind11::module settingsNames = m.def_submodule("opt_settings_names");
  settingsNames.doc() = R"(The ``opt_settings_names`` submodule defines common optimization related setting names to be universal across "
                         "different programs.)";

  pybind11::class_<SettingsNames::Optimizations::Convergence> convergence(settingsNames, "Convergence");
  convergence.def_readonly_static("max_iterations", &SettingsNames::Optimizations::Convergence::maxIter);
  convergence.def_readonly_static("step_max_coefficient", &SettingsNames::Optimizations::Convergence::stepMaxCoeff);
  convergence.def_readonly_static("step_rms", &SettingsNames::Optimizations::Convergence::stepRMS);
  convergence.def_readonly_static("gradient_max_coefficient", &SettingsNames::Optimizations::Convergence::gradMaxCoeff);
  convergence.def_readonly_static("gradient_rms", &SettingsNames::Optimizations::Convergence::gradRMS);
  convergence.def_readonly_static("delta_value", &SettingsNames::Optimizations::Convergence::deltaValue);
  convergence.def_readonly_static("max_iterations", &SettingsNames::Optimizations::Convergence::maxIter);
  convergence.def_readonly_static("requirement", &SettingsNames::Optimizations::Convergence::requirement);

  pybind11::class_<SettingsNames::Optimizations::Bfgs> bfgs(settingsNames, "Bfgs");
  bfgs.def_readonly_static("min_iterations", &SettingsNames::Optimizations::Bfgs::minIter);
  bfgs.def_readonly_static("use_trust_radius", &SettingsNames::Optimizations::Bfgs::useTrustRadius);
  bfgs.def_readonly_static("trust_radius", &SettingsNames::Optimizations::Bfgs::trustRadius);
  bfgs.def_readonly_static("use_gdiis", &SettingsNames::Optimizations::Bfgs::useGdiis);
  bfgs.def_readonly_static("gdiis_max_store", &SettingsNames::Optimizations::Bfgs::gdiisMaxStore);

  pybind11::class_<SettingsNames::Optimizations::Dimer> dimer(settingsNames, "Dimer");
  dimer.def_readonly_static("skip_first_rotation", &SettingsNames::Optimizations::Dimer::skipFirstRotation);
  dimer.def_readonly_static("decrease_rotation_gradient_threshold",
                            &SettingsNames::Optimizations::Dimer::decreaseRotationGradientThreshold);
  dimer.def_readonly_static("gradient_interpolation", &SettingsNames::Optimizations::Dimer::gradientInterpolation);
  dimer.def_readonly_static("rotation_conjugate_gradient", &SettingsNames::Optimizations::Dimer::rotationCG);
  dimer.def_readonly_static("rotation_lbfgs", &SettingsNames::Optimizations::Dimer::rotationLBFGS);
  dimer.def_readonly_static("only_one_rotation", &SettingsNames::Optimizations::Dimer::onlyOneRotation);
  dimer.def_readonly_static("translation", &SettingsNames::Optimizations::Dimer::translation);
  dimer.def_readonly_static("multi_scale", &SettingsNames::Optimizations::Dimer::multiScale);
  dimer.def_readonly_static("radius", &SettingsNames::Optimizations::Dimer::radius);
  dimer.def_readonly_static("phi_tolerance", &SettingsNames::Optimizations::Dimer::phiTolerance);
  dimer.def_readonly_static("rotation_gradient_first", &SettingsNames::Optimizations::Dimer::rotationGradientThresholdFirstCycle);
  dimer.def_readonly_static("rotation_gradient_other",
                            &SettingsNames::Optimizations::Dimer::rotationGradientThresholdOtherCycles);
  dimer.def_readonly_static("lowered_rotation_gradient", &SettingsNames::Optimizations::Dimer::loweredRotationGradientThreshold);
  dimer.def_readonly_static("grad_rmsd_threshold", &SettingsNames::Optimizations::Dimer::gradRMSDthreshold);
  dimer.def_readonly_static("trust_radius", &SettingsNames::Optimizations::Dimer::trustRadius);
  dimer.def_readonly_static("default_translation_step", &SettingsNames::Optimizations::Dimer::defaultTranslationStep);
  dimer.def_readonly_static("max_rotations_first_cycle", &SettingsNames::Optimizations::Dimer::maxRotationsFirstCycle);
  dimer.def_readonly_static("max_rotations_other_cycle", &SettingsNames::Optimizations::Dimer::maxRotationsOtherCycles);
  dimer.def_readonly_static("interval_of_rotations", &SettingsNames::Optimizations::Dimer::intervalOfRotations);
  dimer.def_readonly_static("cycle_of_rotation_gradient_decrease",
                            &SettingsNames::Optimizations::Dimer::cycleOfRotationGradientDecrease);
  dimer.def_readonly_static("lbfgs_memory", &SettingsNames::Optimizations::Dimer::lbfgsMemory);
  dimer.def_readonly_static("bfgs_start", &SettingsNames::Optimizations::Dimer::bfgsStart);
  dimer.def_readonly_static("minimization_cycle", &SettingsNames::Optimizations::Dimer::minimizationCycle);

  pybind11::class_<SettingsNames::Optimizations::Lbfgs> lbfgs(settingsNames, "Lbfgs");
  lbfgs.def_readonly_static("maxm", &SettingsNames::Optimizations::Lbfgs::maxm);
  lbfgs.def_readonly_static("linesearch", &SettingsNames::Optimizations::Lbfgs::linesearch);
  lbfgs.def_readonly_static("c1", &SettingsNames::Optimizations::Lbfgs::c1);
  lbfgs.def_readonly_static("c2", &SettingsNames::Optimizations::Lbfgs::c2);
  lbfgs.def_readonly_static("step_length", &SettingsNames::Optimizations::Lbfgs::stepLength);
  lbfgs.def_readonly_static("use_trust_radius", &SettingsNames::Optimizations::Lbfgs::useTrustRadius);
  lbfgs.def_readonly_static("trust_radius", &SettingsNames::Optimizations::Lbfgs::trustRadius);
  lbfgs.def_readonly_static("max_backtracking", &SettingsNames::Optimizations::Lbfgs::maxBacktracking);

  pybind11::class_<SettingsNames::Optimizations::SteepestDescent> steepestdescent(settingsNames, "SteepestDescent");
  steepestdescent.def_readonly_static("factor", &SettingsNames::Optimizations::SteepestDescent::factor);
  steepestdescent.def_readonly_static("use_trust_radius", &SettingsNames::Optimizations::SteepestDescent::useTrustRadius);
  steepestdescent.def_readonly_static("trust_radius", &SettingsNames::Optimizations::SteepestDescent::trustRadius);
  steepestdescent.def_readonly_static("dynamic_multiplier", &SettingsNames::Optimizations::SteepestDescent::dynamicMultiplier);

  pybind11::class_<SettingsNames::Optimizations::Bofill> bofill(settingsNames, "Bofill");
  bofill.def_readonly_static("trust_radius", &SettingsNames::Optimizations::Bofill::trustRadius);
  bofill.def_readonly_static("hessian_update", &SettingsNames::Optimizations::Bofill::hessianUpdate);
  bofill.def_readonly_static("follow_mode", &SettingsNames::Optimizations::Bofill::mode);

  pybind11::class_<SettingsNames::Optimizations::EigenvectorFollowing> eigenvectorfollowing(settingsNames,
                                                                                            "EigenvectorFollowing");
  eigenvectorfollowing.def_readonly_static("trust_radius", &SettingsNames::Optimizations::EigenvectorFollowing::trustRadius);
  eigenvectorfollowing.def_readonly_static("follow_mode", &SettingsNames::Optimizations::EigenvectorFollowing::mode);

  pybind11::class_<SettingsNames::Optimizations::NewtonRaphson> newtonraphson(settingsNames, "NewtonRaphson");
  newtonraphson.def_readonly_static("trust_radius", &SettingsNames::Optimizations::NewtonRaphson::trustRadius);
  newtonraphson.def_readonly_static("svd_threshold", &SettingsNames::Optimizations::NewtonRaphson::svdThreshold);

  pybind11::class_<SettingsNames::Optimizations::GeometryOptimizer> geometryOptimizer(settingsNames,
                                                                                      "GeometryOptimizer");
  geometryOptimizer.def_readonly_static("fixed_atoms", &SettingsNames::Optimizations::GeometryOptimizer::fixedAtoms);
  geometryOptimizer.def_readonly_static("coordinate_system", &SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem);

  pybind11::class_<SettingsNames::Optimizations::CellOptimizer> cellOptimizer(settingsNames, "CellOptimizer");
  cellOptimizer.def_readonly_static("optimize_angles", &SettingsNames::Optimizations::CellOptimizer::optimizeAngles);
  cellOptimizer.def_readonly_static("optimize_a", &SettingsNames::Optimizations::CellOptimizer::optimizeA);
  cellOptimizer.def_readonly_static("optimize_b", &SettingsNames::Optimizations::CellOptimizer::optimizeB);
  cellOptimizer.def_readonly_static("optimize_c", &SettingsNames::Optimizations::CellOptimizer::optimizeC);
  cellOptimizer.def_readonly_static("geoopt_max_convergence_iterations",
                                    &SettingsNames::Optimizations::CellOptimizer::geooptMaxIterations);
  cellOptimizer.def_readonly_static("cellopt_max_convergence_iterations",
                                    &SettingsNames::Optimizations::CellOptimizer::celloptMaxIterations);

  pybind11::class_<SettingsNames::Optimizations::Nt> nt(settingsNames, "Nt");
  nt.def_readonly_static("rhs_list", &SettingsNames::Optimizations::Nt::rHSList);
  nt.def_readonly_static("lhs_list", &SettingsNames::Optimizations::Nt::lHSList);
  nt.def_readonly_static("attractive", &SettingsNames::Optimizations::Nt::attractive);
  nt.def_readonly_static("total_force_norm", &SettingsNames::Optimizations::Nt::totalForceNorm);
  nt.def_readonly_static("convergence_max_iterations", &SettingsNames::Optimizations::Nt::maxIter);
  nt.def_readonly_static("convergence_repulsive_stop", &SettingsNames::Optimizations::Nt::repulsiveStop);
  nt.def_readonly_static("convergence_attractive_stop", &SettingsNames::Optimizations::Nt::attractiveStop);
  nt.def_readonly_static("sd_factor", &SettingsNames::Optimizations::Nt::sdFactor);
  nt.def_readonly_static("use_micro_cycles", &SettingsNames::Optimizations::Nt::useMicroCycles);
  nt.def_readonly_static("fixed_number_of_micro_cycles", &SettingsNames::Optimizations::Nt::fixedNumberOfMicroCycles);
  nt.def_readonly_static("number_of_micro_cycles", &SettingsNames::Optimizations::Nt::numberOfMicroCycles);
  nt.def_readonly_static("filter_passes", &SettingsNames::Optimizations::Nt::filterPasses);
  nt.def_readonly_static("extraction_criterion", &SettingsNames::Optimizations::Nt::extractionCriterion);
  nt.def_readonly_static("coordinate_system", &SettingsNames::Optimizations::Nt::coordinateSystem);
  nt.def_readonly_static("constrained_atoms", &SettingsNames::Optimizations::Nt::fixedAtoms);
  nt.def_readonly_static("movable_side", &SettingsNames::Optimizations::Nt::movableSide);
  nt.def_readonly_static("highest_maximum", &SettingsNames::Optimizations::Nt::extractHighest);
  nt.def_readonly_static("first_maximum", &SettingsNames::Optimizations::Nt::extractFirst);

  pybind11::class_<SettingsNames::Optimizations::Nt2> nt2(settingsNames, "Nt2");
  nt2.def_readonly_static("associations", &SettingsNames::Optimizations::Nt2::assList);
  nt2.def_readonly_static("dissociations", &SettingsNames::Optimizations::Nt2::dissList);
  nt2.def_readonly_static("total_force_norm", &SettingsNames::Optimizations::Nt2::totalForceNorm);
  nt2.def_readonly_static("convergence_max_iterations", &SettingsNames::Optimizations::Nt2::maxIter);
  nt2.def_readonly_static("convergence_attractive_stop", &SettingsNames::Optimizations::Nt2::attractiveStop);
  nt2.def_readonly_static("sd_factor", &SettingsNames::Optimizations::Nt2::sdFactor);
  nt2.def_readonly_static("use_micro_cycles", &SettingsNames::Optimizations::Nt2::useMicroCycles);
  nt2.def_readonly_static("fixed_number_of_micro_cycles", &SettingsNames::Optimizations::Nt2::fixedNumberOfMicroCycles);
  nt2.def_readonly_static("number_of_micro_cycles", &SettingsNames::Optimizations::Nt2::numberOfMicroCycles);
  nt2.def_readonly_static("filter_passes", &SettingsNames::Optimizations::Nt2::filterPasses);
  nt2.def_readonly_static("extraction_criterion", &SettingsNames::Optimizations::Nt2::extractionCriterion);
  nt2.def_readonly_static("coordinate_system", &SettingsNames::Optimizations::Nt2::coordinateSystem);
  nt2.def_readonly_static("constrained_atoms", &SettingsNames::Optimizations::Nt2::fixedAtoms);
  nt2.def_readonly_static("last_maximum_before_first_target", &SettingsNames::Optimizations::Nt2::extractLastBeforeTarget);
  nt2.def_readonly_static("highest_maximum", &SettingsNames::Optimizations::Nt2::extractHighest);
  nt2.def_readonly_static("first_maximum", &SettingsNames::Optimizations::Nt2::extractFirst);

  pybind11::class_<SettingsNames::Optimizations::Irc> irc(settingsNames, "Irc");
  irc.def_readonly_static("initial_step_size", &SettingsNames::Optimizations::Irc::initialStepSize);
  irc.def_readonly_static("coordinate_system", &SettingsNames::Optimizations::Irc::coordinateSystem);

  pybind11::class_<SettingsNames::Optimizations::Afir> afir(settingsNames, "Afir");
  afir.def_readonly_static("rhs_list", &SettingsNames::Optimizations::Afir::rHSList);
  afir.def_readonly_static("lhs_list", &SettingsNames::Optimizations::Afir::lHSList);
  afir.def_readonly_static("weak_forces", &SettingsNames::Optimizations::Afir::weakForces);
  afir.def_readonly_static("attractive", &SettingsNames::Optimizations::Afir::attractive);
  afir.def_readonly_static("energy_allowance", &SettingsNames::Optimizations::Afir::energyAllowance);
  afir.def_readonly_static("phase_in", &SettingsNames::Optimizations::Afir::phaseIn);
  afir.def_readonly_static("coordinate_system", &SettingsNames::Optimizations::Afir::coordinateSystem);
}
