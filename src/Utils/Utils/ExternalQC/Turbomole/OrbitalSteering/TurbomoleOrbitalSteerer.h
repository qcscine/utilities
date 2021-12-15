/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef ORBITALPERTURBATION_TURBOMOLESTEERER_H
#define ORBITALPERTURBATION_TURBOMOLESTEERER_H

#include <Utils/ExternalQC/Turbomole/TurbomoleFiles.h>
#include <string>
#include <vector>

namespace Scine {
namespace Utils {
class MolecularOrbitals;
namespace ExternalQC {
namespace TurbomoleOrbitalPerturbation {

class TurbomoleOrbitalSteerer {
 public:
  /**
   * @brief Constructor
   */
  TurbomoleOrbitalSteerer(std::string& calculationDirectory);
  /**
   * @brief Performs the orbital steering.
   */
  void steerOrbitals();

 private:
  /*
   * @brief Parses number of alpha and beta electrons from the control file.
   */
  std::pair<int, int> extractNumberOfAlphaAndBetaElectrons();
  /*
   * @brief Parses the total number of orbitals from the control file.
   */
  int extractNumerOfOrbitals();
  /*
   * @brief Mixes the orbitals.
   */
  void mixOrbitals(MolecularOrbitals& mo, unsigned nAlpha, unsigned nBeta);

  std::string& calculationDirectory_;
  TurbomoleFiles files_;
};

} // namespace TurbomoleOrbitalPerturbation
} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // ORBITALPERTURBATION_TURBOMOLESTEERER_H