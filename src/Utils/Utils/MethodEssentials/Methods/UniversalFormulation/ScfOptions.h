/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SCFOPTIONS_H
#define UTILS_SCFOPTIONS_H

#include <Utils/MethodEssentials/MethodFactories/MixerFactory.h>
#include <string>

namespace Scine {

namespace Utils {
namespace UniversalSettings {
class DescriptorCollection;
class ValueCollection;
} // namespace UniversalSettings
} // namespace Utils

namespace Utils {

class SCFMethod;
enum class scf_mixer_t;

/*!
 * SCF options for the UniversalSettings syntax.
 */
struct ScfOptions {
  static constexpr const char* name = "scf_options";
  static constexpr const char* description = "SCF settings";
  static constexpr const char* mixerName = "scf_mixer";
  static constexpr const char* mixerDescription = "SCF convergence accelerator";
  static constexpr const char* thresholdName = "scf_threshold";
  static constexpr const char* thresholdDescription = "SCF threshold (RMSD of density matrix)";
  static constexpr const char* maxIterationName = "scf_max_iterations";
  static constexpr const char* maxIterationDescription = "Maximal number of SCF iterations";

  // Mixer options:
  static constexpr const char* noMixer = "none";
  static constexpr const char* diis = "DIIS";
  static constexpr const char* ediis = "EDIIS";
  static constexpr const char* ediisDiis = "EDIIS+DIIS";

  static Scine::Utils::UniversalSettings::DescriptorCollection getSettingDescriptors();
  static void applySettings(const Scine::Utils::UniversalSettings::ValueCollection& values, SCFMethod& method);
  static Scine::Utils::UniversalSettings::ValueCollection getAppliedSettings(const SCFMethod& method);

 private:
  static scf_mixer_t convertToMixer(const std::string& s);
  static std::string convertToString(scf_mixer_t m);
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_SCFOPTIONS_H