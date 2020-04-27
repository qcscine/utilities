/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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
   * @brief Get the covalent radius for a given element.
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
