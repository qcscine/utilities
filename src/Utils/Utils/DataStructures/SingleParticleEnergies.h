/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SINGLEPARTICLEENERGIES_H
#define UTILS_SINGLEPARTICLEENERGIES_H

#include <Eigen/Core>
#include <cassert>
#include <vector>

namespace Scine {
namespace Utils {

/*!
 * Class handling the single-particle energies obtained as
 * eigenvalues of the (generalized) eigenvalue problem.
 * NB: no check is made if the correct (i.e., Restricted / Unrestricted)
 * function variants are called.
 */

class SingleParticleEnergies {
 public:
  using EnergyLevels = std::vector<double>;

  static SingleParticleEnergies createEmptyUnrestrictedEnergies();
  static SingleParticleEnergies createEmptyRestrictedEnergies();

  bool isRestricted() const;

  int getRestrictedNLevels() const;
  int getUnrestrictedNLevels() const;

  void setRestricted(const Eigen::VectorXd& values);
  void setUnrestricted(const Eigen::VectorXd& alpha, const Eigen::VectorXd& beta);

  double getRestrictedLevelEnergy(int index) const;
  double getAlphaLevelEnergy(int index) const;
  double getBetaLevelEnergy(int index) const;

  const EnergyLevels& getRestrictedEnergies() const;
  const EnergyLevels& getAlphaEnergies() const;
  const EnergyLevels& getBetaEnergies() const;

  bool restricted_{true};
  EnergyLevels restrictedEnergies_;
  EnergyLevels alphaEnergies_;
  EnergyLevels betaEnergies_;
};

inline const SingleParticleEnergies::EnergyLevels& SingleParticleEnergies::getBetaEnergies() const {
  return betaEnergies_;
}

inline const SingleParticleEnergies::EnergyLevels& SingleParticleEnergies::getAlphaEnergies() const {
  return alphaEnergies_;
}

inline const SingleParticleEnergies::EnergyLevels& SingleParticleEnergies::getRestrictedEnergies() const {
  return restrictedEnergies_;
}

inline double SingleParticleEnergies::getRestrictedLevelEnergy(int index) const {
  assert(restricted_ && "Cannot return restricted level energy because calculation was unrestricted.");
  assert(0 <= index && index < static_cast<int>(restrictedEnergies_.size()));
  return restrictedEnergies_[index];
}

inline double SingleParticleEnergies::getAlphaLevelEnergy(int index) const {
  assert(!restricted_ && "Cannot return unrestricted level energy because calculation was restricted.");
  assert(0 <= index && index < static_cast<int>(alphaEnergies_.size()));
  return alphaEnergies_[index];
}

inline double SingleParticleEnergies::getBetaLevelEnergy(int index) const {
  assert(!restricted_ && "Cannot return unrestricted level energy because calculation was restricted.");
  assert(0 <= index && index < static_cast<int>(betaEnergies_.size()));
  return betaEnergies_[index];
}

} // namespace Utils
} // namespace Scine
#endif // UTILS_SINGLEPARTICLEENERGIES_H
