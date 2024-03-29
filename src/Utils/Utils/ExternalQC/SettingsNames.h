/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_SETTINGSNAMES_H
#define UTILS_EXTERNALQC_SETTINGSNAMES_H

namespace Scine {
namespace Utils {
namespace ExternalQC {
namespace SettingsNames {

// global
static constexpr const char* baseWorkingDirectory = "base_working_directory";
static constexpr const char* deleteTemporaryFiles = "delete_tmp_files";
static constexpr const char* steerOrbitals = "steer_orbitals";
static constexpr const char* pointChargesFile = "point_charges_file";
static constexpr const char* enforceScfCriterion = "enforce_scf_criterion";
static constexpr const char* hessianCalculationType = "hessian_calculation_type";
// orca
static constexpr const char* orcaFilenameBase = "orca_filename_base";
static constexpr const char* specialOption = "special_option";
static constexpr const char* performBrokenSymmetryCalculation = "perform_broken_symmetry_calculation";
static constexpr const char* initialSpinMultiplicity = "initial_spin_multiplicity";
static constexpr const char* spinFlipSites = "spin_flip_sites";
static constexpr const char* calculateMoessbauerParameter = "calculate_moessbauer";
static constexpr const char* orcaAuxCBasisSet = "auxc_basis_set";
static constexpr const char* orcaCabsBasisSet = "cabs_basis_set";
static constexpr const char* gradientCalculationType = "gradient_calculation_type";
// cp2k
static constexpr const char* cp2kFilenameBase = "cp2k_filename_base";
static constexpr const char* moleculePeriodicBoundaries = "28.35, 28.35, 28.35, 90.0, 90.0, 90.0, XYZ";
static constexpr const char* relMultiGridCutoff = "relative_multi_grid_cutoff";
static constexpr const char* orbitalTransformation = "orbital_transformation";
static constexpr const char* outerScf = "outer_scf";
static constexpr const char* poissonSolver = "poisson_solver";
static constexpr const char* additionalMos = "additional_mos";
static constexpr const char* allowUnconvergedScf = "allow_unconverged_scf";
static constexpr const char* planeWaveCutoff = "plane_wave_cutoff";
static constexpr const char* nGrids = "n_grids";
static constexpr const char* scfGuess = "scf_guess";
static constexpr const char* dipoleCorrection = "dipole_correction";
static constexpr const char* additionalOutputFile = "additional_output_file";
// gaussian
static constexpr const char* gaussianFilenameBase = "gaussian_filename_base";
// turbomole
static constexpr const char* scfOrbitalShift = "scf_orbitalshift";
static constexpr const char* scfDampingValue = "scf_damping_value";
static constexpr const char* enableRi = "enable_ri";
static constexpr const char* numExcitedStates = "num_excited_states";
static constexpr const char* dftGrid = "dft_grid";
static constexpr const char* cavityPointsPerAtom = "cavity_points_per_atom";
static constexpr const char* cavitySegmentsPerAtom = "cavity_segments_per_atom";
static constexpr const char* enforceNumforce = "enforce_numforce";

} // namespace SettingsNames
} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_SETTINGSNAMES_H
