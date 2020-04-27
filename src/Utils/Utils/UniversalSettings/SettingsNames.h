/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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
static constexpr const char* unrestrictedCalculation = "unrestricted_calculation";
static constexpr const char* selfConsistanceCriterion = "self_consistence_criterion";
static constexpr const char* maxIterations = "max_scf_iterations";
static constexpr const char* mixer = "scf_mixer";
static constexpr const char* loggerVerbosity = "log";
static constexpr const char* temperature = "temperature";
static constexpr const char* davidsonForGroundState = "davidson_for_ground_state";

static constexpr const char* NDDODipoleApproximation = "nddo_dipole";
static constexpr const char* parameterFile = "parameter_file";
static constexpr const char* parameterRootDirectory = "parameter_root";

static constexpr const char* method = "method";
static constexpr const char* basisSet = "basis_set";

//! @brief Struct to contain the name of the mixers available.
struct ScfMixers {
  static constexpr const char* noMixer = "no_mixer";
  static constexpr const char* diis = "diis";
  static constexpr const char* ediis = "ediis";
  static constexpr const char* ediisDiis = "ediis_diis";
};

//! @brief Struct to contain the name of the level of verbosity for the logger.
struct LogLevels {
  static constexpr const char* none = "none";
  static constexpr const char* trace = "trace";
  static constexpr const char* debug = "debug";
  static constexpr const char* info = "info";
  static constexpr const char* warning = "warning";
  static constexpr const char* error = "error";
  static constexpr const char* fatal = "fatal";
};

//! @@brief Settings for linear response time dependent methods.
static constexpr const char* numberOfEigenstates = "number_eigenstates";
static constexpr const char* initialSubspaceDimension = "initial_subspace_dimension";
static constexpr const char* useSparseImplementation = "LRTD_sparse_implementation";
static constexpr const char* spinBlock = "spin_block";
struct SpinBlocks {
  static constexpr const char* singlet = "singlet";
  static constexpr const char* triplet = "triplet";
  static constexpr const char* singletAndTriplet = "both";
};

static constexpr const char* directness = "directness";
struct Directness {
  static constexpr const char* direct = "direct";
  static constexpr const char* standard = "standard";
};

} // namespace SettingsNames
} // namespace Utils
} // namespace Scine

#endif // UTILS_SETTINGSNAMES_H
