/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_AUFBAUPRINCIPLEOCCUPATIONGENERATOR_H
#define UTILS_AUFBAUPRINCIPLEOCCUPATIONGENERATOR_H

#include <Utils/MethodEssentials/util/LcaoUtil/ElectronicOccupationGenerator.h>

namespace Scine {
namespace Utils {

namespace LcaoUtil {

class AufbauPrincipleOccupationGenerator : public ElectronicOccupationGenerator {
 public:
 private:
  ElectronicOccupation generateOccupationImpl() override;
};

} // namespace LcaoUtil

} // namespace Utils
} // namespace Scine
#endif // UTILS_AUFBAUPRINCIPLEOCCUPATIONGENERATOR_H
