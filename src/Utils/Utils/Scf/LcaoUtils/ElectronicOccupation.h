/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_ELECTRONICOCCUPATION_H
#define UTILS_ELECTRONICOCCUPATION_H

#include <vector>

namespace Scine {
namespace Utils {
namespace LcaoUtils {

/*!
 * Class to hold information about which molecular orbitals are occupied.
 * TODO: Implement fractional occupation ?
 */
class ElectronicOccupation {
 public:
  int numberRestrictedElectrons() const;
  int numberOccupiedRestrictedOrbitals() const;
  int numberAlphaElectrons() const;
  int numberBetaElectrons() const;
  /*! Returns true if the electrons all occupy the lowest orbitals. */
  bool isFilledUpFromTheBottom() const;
  bool hasUnpairedRHFElectron() const;

  bool isRestricted() const;
  bool isUnrestricted() const;
  /*! From a restricted occupation, transform to an unrestricted occupation. */
  void makeUnrestricted();
  ElectronicOccupation toUnrestricted() const;

  const std::vector<int>& getFilledRestrictedOrbitals() const;
  const std::vector<int>& getFilledAlphaOrbitals() const;
  const std::vector<int>& getFilledBetaOrbitals() const;

  void fillLowestRestrictedOrbitalsWithElectrons(int nElectrons);
  void fillLowestUnrestrictedOrbitals(int nAlphaElectrons, int nBetaElectrons);
  void fillSpecifiedRestrictedOrbitals(std::vector<int> occupiedRestrictedOrbitals);
  void fillSpecifiedUnrestrictedOrbitals(std::vector<int> occupiedAlphaOrbitals, std::vector<int> occupiedBetaOrbitals);

 private:
  void reset();
  void checkWhetherFilledFromBottom() const;
  std::vector<int> generateOccupiedArrayFromNumberOfOccupiedOrbitals(int nOrbitals) const;

  bool restricted_ = false;
  int nRestrictedElectrons_ = 0;
  int nAlphaElectrons_ = 0;
  int nBetaElectrons_ = 0;

  bool hasUnpairedElectron_ = false;
  mutable bool specifiedWithNumberElectrons_ = false;
  mutable bool specifiedWithOrbitals_ = false;

  mutable std::vector<int> filledRestrictedOrbitals_;
  mutable std::vector<int> filledAlphaOrbitals_;
  mutable std::vector<int> filledBetaOrbitals_;
};

} // namespace LcaoUtils
} // namespace Utils
} // namespace Scine

#endif // UTILS_ELECTRONICOCCUPATION_H
