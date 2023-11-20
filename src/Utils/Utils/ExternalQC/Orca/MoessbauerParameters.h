/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ORCAMOESSBAUERCONTAINER_H
#define UTILS_ORCAMOESSBAUERCONTAINER_H

#include "Utils/ExternalQC/SettingsNames.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Settings.h"
#include <Utils/Geometry/ElementInfo.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {

namespace Moessbauer {
struct MoessbauerParameterContainer {
  int numIrons = 0;
  std::vector<double> quadrupoleSplittings;
  std::vector<double> etas;
  std::vector<double> densities;

  bool isApprox(const MoessbauerParameterContainer& rhs) const {
    return (std::fabs(this->numIrons - rhs.numIrons) == 0 &&
            std::equal(this->quadrupoleSplittings.begin(), this->quadrupoleSplittings.end(), rhs.quadrupoleSplittings.begin()) &&
            std::equal(this->etas.begin(), this->etas.end(), rhs.etas.begin()) &&
            std::equal(this->densities.begin(), this->densities.end(), rhs.densities.begin()));
  }
};

inline bool moessbauerNeededAndPossible(const Utils::AtomCollection& structure, const Settings& settings) {
  bool moessbauerRequired = settings.getBool(Utils::ExternalQC::SettingsNames::calculateMoessbauerParameter);
  bool structureContainsIron = std::any_of(std::begin(structure), std::end(structure), [&](const Utils::Atom& a) {
    return (a.getElementType() == Utils::ElementType::Fe);
  });
  return (moessbauerRequired && structureContainsIron);
}

inline int determineNumIrons(const Utils::AtomCollection& structure) {
  return std::count_if(structure.begin(), structure.end(),
                       [&](const Utils::Atom& a) { return a.getElementType() == Utils::ElementType::Fe; });
}

} // namespace Moessbauer
} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_ORCACALCULATOR_H
