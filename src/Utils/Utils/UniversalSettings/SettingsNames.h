/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SETTINGSNAMES_H
#define UTILS_SETTINGSNAMES_H

namespace Scine {
namespace Utils {
namespace SettingsNames {

//! SettingsNames provides a consistent list of settings names throughout the whole module.

static constexpr const char* molecularCharge = "molecular_charge";
static constexpr const char* spinMultiplicity = "spin_multiplicity";
static constexpr const char* maxScfIterations = "max_scf_iterations";
static constexpr const char* selfConsistenceCriterion = "self_consistence_criterion";
static constexpr const char* densityRmsdCriterion = "density_rmsd_criterion";
static constexpr const char* mixer = "scf_mixer";
static constexpr const char* loggerVerbosity = "log";
static constexpr const char* symmetryNumber = "symmetry_number";
static constexpr const char* methodParameters = "method_parameters";
static constexpr const char* NDDODipoleApproximation = "nddo_dipole";
static constexpr const char* mmCharges = "mm_charges";
static constexpr const char* parameterFilePath = "mm_parameter_file";

// Model
static constexpr const char* method = "method";
static constexpr const char* methodFamily = "method_family";
static constexpr const char* spinMode = "spin_mode";
static constexpr const char* program = "program";
static constexpr const char* version = "version";
static constexpr const char* basisSet = "basis_set";
static constexpr const char* temperature = "temperature";
static constexpr const char* pressure = "pressure";
static constexpr const char* electronicTemperature = "electronic_temperature";
static constexpr const char* solvation = "solvation";
static constexpr const char* solvent = "solvent";
static constexpr const char* embedding = "embedding";
static constexpr const char* periodicBoundaries = "periodic_boundaries";
static constexpr const char* externalField = "external_field";

static constexpr const char* externalProgramMemory = "external_program_memory";
static constexpr const char* externalProgramNProcs = "external_program_nprocs";

static constexpr const char* scfDamping = "scf_damping";

// QM/MM
static constexpr const char* qmAtomsList = "qm_atoms";
static constexpr const char* electrostaticEmbedding = "electrostatic_embedding";
static constexpr const char* ignoreQmOption = "ignore_qm";
static constexpr const char* optimizeLinks = "optimize_links";

//! @brief Struct to contain the name of the mixers available.
struct ScfMixers {
  static constexpr const char* noMixer = "no_mixer";
  static constexpr const char* diis = "diis";
  static constexpr const char* ediis = "ediis";
  static constexpr const char* ediisDiis = "ediis_diis";
};

//! @@brief Settings for linear response time dependent methods.
static constexpr const char* maxDavidsonIterations = "max_davidson_iterations";
//! @@brief Path to the excited states parameter file
static constexpr const char* excitedStatesParamFile = "excited_parameterfile";
//! @brief The number of roots to be calculated.
static constexpr const char* numberOfEigenstates = "number_eigenstates";
//! @brief The initial guess space, bigger spaces could speed up convergence.
static constexpr const char* initialSubspaceDimension = "initial_subspace_dimension";
//! @brief The maximal memory allowed by the program.
static constexpr const char* maxMemory = "max_memory";
//! @brief Defines the spin block to calculate in excited states calculation.
static constexpr const char* spinBlock = "spin_block";
struct SpinBlocks {
  static constexpr const char* singlet = "singlet";
  static constexpr const char* triplet = "triplet";
  static constexpr const char* singletAndTriplet = "both";
};

//! @brief Sets whether the basis of singly excited determinants should be pruned and with which method.
static constexpr const char* pruneBasis = "prune_basis";
struct PruningOptions {
  static constexpr const char* none = "none";
  static constexpr const char* energy = "energy";
};
//! @brief Sets the threshold for pruning with an energy criterion in au.
static constexpr const char* energyThreshold = "energy_threshold";
//! @brief Sets the threshold for pruning with an intensity criterion in au.
static constexpr const char* perturbativeThreshold = "pt_threshold";

//! @brief Davidson options
static constexpr const char* directness = "directness";
struct Directness {
  static constexpr const char* direct = "direct";
  static constexpr const char* standard = "standard";
};

} // namespace SettingsNames
} // namespace Utils
} // namespace Scine

#endif // UTILS_SETTINGSNAMES_H
