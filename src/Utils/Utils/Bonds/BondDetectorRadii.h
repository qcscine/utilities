/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_BONDDETECTORRADII_H_
#define UTILS_BONDDETECTORRADII_H_

#include <Utils/Constants.h>
#include <Utils/Geometry/ElementTypes.h>
#include <array>

namespace Scine {
namespace Utils {

/**
 * @class BondDetectorRadii BondDetectorRadii.h
 * @brief Holds the covalent radii required in BondDetector.
 */
class BondDetectorRadii {
 public:
  /// @brief Default constructor.
  BondDetectorRadii();
  /**
   * @brief Get the Cambridge Structural Database covalent radius for a given element.
   *
   * For molecular structures the covalent radii available in Utils::ElementInfo are often more appropriate.
   *
   * The covalent radii here were extracted from the Cambridge Structural Database (CSD) on 04/08/2020:\n
   * https://www.ccdc.cam.ac.uk/support-and-resources/ccdcresources/Elemental_Radii.xlsx\n
   * C. R. Groom, I. J. Bruno, M. P. Lightfoot and S. C. Ward, Acta Cryst. 2016, B72, 171-179.
   * [DOI 10.1107/S2052520616003954]\n
   *
   * The original reference for the covalent radii of metals is:\n
   * S. Alvarez et al., Dalton Trans., 2008, 2832-2838. [DOI: 10.1039/B801115J]
   *
   * Elements not yet encountered in the CSD get assigned a default covalent radius of 1.5 Angstrom.
   *
   * @param e The element.
   * @return double The covalent radius.
   */
  double getRadius(ElementType e) const;

 private:
  void fillArray();
  void setRadiusInAngstrom(ElementType e, double rInAngstrom);
  void setRadius(ElementType e, double r);

  static constexpr int numberElements = 120;
  static constexpr double defaultRadius = toBohr(Angstrom(1.50));
  std::array<double, numberElements> radii;
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_BONDDETECTORRADII_H_
