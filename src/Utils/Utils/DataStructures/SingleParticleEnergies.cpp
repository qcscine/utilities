/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SingleParticleEnergies.h"

namespace Scine {
namespace Utils {

namespace {
inline void assignStdVectorFromEigenVectorXd(std::vector<double>& v1, const Eigen::VectorXd& v2) {
  v1.assign(v2.data(), v2.data() + v2.size());
}
} // namespace

bool SingleParticleEnergies::isRestricted() const {
  return restricted_;
}

int SingleParticleEnergies::getRestrictedNLevels() const {
  return static_cast<int>(restrictedEnergies_.size());
}

int SingleParticleEnergies::getUnrestrictedNLevels() const {
  return static_cast<int>(alphaEnergies_.size());
}

void SingleParticleEnergies::setRestricted(const Eigen::VectorXd& values) {
  restricted_ = true;
  alphaEnergies_.clear();
  betaEnergies_.clear();
  assignStdVectorFromEigenVectorXd(restrictedEnergies_, values);
}

void SingleParticleEnergies::setUnrestricted(const Eigen::VectorXd& alpha, const Eigen::VectorXd& beta) {
  restricted_ = false;
  restrictedEnergies_.clear();
  assignStdVectorFromEigenVectorXd(alphaEnergies_, alpha);
  assignStdVectorFromEigenVectorXd(betaEnergies_, beta);
}

SingleParticleEnergies SingleParticleEnergies::createEmptyUnrestrictedEnergies() {
  SingleParticleEnergies energies;
  energies.setUnrestricted(Eigen::VectorXd{}, Eigen::VectorXd{});
  return energies;
}

SingleParticleEnergies SingleParticleEnergies::createEmptyRestrictedEnergies() {
  SingleParticleEnergies energies;
  energies.setRestricted(Eigen::VectorXd{});
  return energies;
}

} // namespace Utils
} // namespace Scine
