/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SCFCONVERGENCEACCELERATOR_H
#define UTILS_SCFCONVERGENCEACCELERATOR_H

#include <memory>

namespace Scine {
namespace Utils {

class ScfMethod;
class ScfModifier;
enum class scf_mixer_t;

/*!
 * This class sets up the convergence acceleration for a SCF method.
 */
class ScfConvergenceAccelerator {
 public:
  explicit ScfConvergenceAccelerator(ScfMethod& method);
  ~ScfConvergenceAccelerator();

  void setScfMixer(scf_mixer_t mixer);
  scf_mixer_t getScfMixer() const;

 private:
  void removeCurrentMixer();
  void setMixer(scf_mixer_t mixer);

  ScfMethod& method_;
  scf_mixer_t currentScheme_;
  std::shared_ptr<ScfModifier> activeScfMixer_;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_SCFCONVERGENCEACCELERATOR_H