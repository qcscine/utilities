/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef ORBITALPERTURBATION_TURBOMOLEORBITALSMETAINFORMATION_H
#define ORBITALPERTURBATION_TURBOMOLEORBITALSMETAINFORMATION_H

#include <string>
#include <vector>

namespace Scine {
namespace Utils {
namespace ExternalQC {
namespace TurbomoleOrbitalPerturbation {

/*!
 * Information contained in a turbomole molecular orbitals file apart from the molecular orbitals themselves.
 */
class TurbomoleOrbitalsMetaInformation {
 public:
  void setHeader(std::string h1, std::string h2, std::string h3);
  void setFooter(std::string footer);
  void addOrbitalInformation(std::string orbitalInformation);

  const std::string& getFirstHeaderLine() const;
  const std::string& getSecondHeaderLine() const;
  const std::string& getThirdHeaderLine() const;
  const std::string& getFooter() const;

  const std::string& getOrbitalInformation(unsigned index) const;

 private:
  std::string header1_, header2_, header3_;
  std::string footer_;
  std::vector<std::string> orbitalInformation_;
};

} // namespace TurbomoleOrbitalPerturbation
} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // ORBITALPERTURBATION_TURBOMOLEORBITALSMETAINFORMATION_H