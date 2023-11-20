/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_COORDINATESYSTEM_H_
#define UTILS_COORDINATESYSTEM_H_

#include <stdexcept>
#include <string>

namespace Scine {
namespace Utils {

/**
 * @brief Enum Class specifying the applied coordinate system
 *
 * Internal represents true internal coordinates defined by bond length, angle, dihedral etc.
 * CartesianWithoutRotTrans represents Cartesian coordinates with the degrees of freedom in translation and rotation
 * removed
 * Cartesian represents true Cartesian coordinates with 3N parameters for N atoms
 */
enum class CoordinateSystem {
  Internal,                 // true internal coordinates
  CartesianWithoutRotTrans, // Cartesian without translation / rotation
  Cartesian                 // true Cartesian coordinates
};

class CoordinateSystemInterpreter {
 public:
  // @brief static functions only
  CoordinateSystemInterpreter() = delete;
  static CoordinateSystem getCoordinateSystemFromString(const std::string& coordinateSystem) {
    if (coordinateSystem == "internal") {
      return CoordinateSystem::Internal;
    }
    if (coordinateSystem == "cartesianWithoutRotTrans") {
      return CoordinateSystem::CartesianWithoutRotTrans;
    }
    if (coordinateSystem == "cartesian") {
      return CoordinateSystem::Cartesian;
    }
    throw std::logic_error("Unknown coordinate system " + coordinateSystem);
  };

  static std::string getStringFromCoordinateSystem(const CoordinateSystem& coordinateSystem) {
    if (coordinateSystem == CoordinateSystem::Internal) {
      return "internal";
    }
    if (coordinateSystem == CoordinateSystem::CartesianWithoutRotTrans) {
      return "cartesianWithoutRotTrans";
    }
    if (coordinateSystem == CoordinateSystem::Cartesian) {
      return "cartesian";
    }
    throw std::logic_error("Unknown string representation for this coordinate system.");
  };
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_COORDINATESYSTEM_H_
