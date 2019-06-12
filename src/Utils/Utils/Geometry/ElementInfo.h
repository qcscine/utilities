/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ELEMENTINFO_H
#define UTILS_ELEMENTINFO_H

#include "Utils/Geometry/ElementTypes.h"
#include <map>
#include <stdexcept>
#include <string>

namespace Scine {
namespace Utils {

/**
 * @class ElementSymbolNotFound ElementInfo.h
 * @brief A runtime error specific to an element not being found.
 */
class ElementSymbolNotFound : public std::runtime_error {
 public:
  explicit ElementSymbolNotFound(const std::string& symbol) : runtime_error(symbol + " is not a known element.") {
  }
};

/**
 * @class ElementInfo ElementInfo.h
 * @brief Provides information about elements, such as mass, van-der-Waals radius, etc.
 *
 * This class only wraps around the actual data and their handling.
 * For the underlying data see ElementInfoData.h and ElementInfoData.cpp.
 */
class ElementInfo {
 public:
  /**
   * @brief Translation from std::string to ElementType enum.
   * @param symbol The atom symbol. First character captitalized, other lower case.
   * @return ElementType Returns the ElementType.
   */
  static ElementType elementTypeForSymbol(const std::string& symbol);
  /**
   * @brief Translation from ElementType enum to std::string.
   * @param e The ElementType.
   * @return std::string Returns the symbol string.
   */
  static std::string symbol(ElementType e);
  /**
   * @brief Getter for the mass (isotope average, precision: 3 digits).
   * @param e The ElementType.
   * @return double Returns the mass of the element.
   */
  static double mass(ElementType e);
  /**
   * @brief Getter for the van der Waals radius in atomic units.
   * @param e The ElementType.
   * @return double Returns the van der Waals radius in atomic units.
   */
  static double vdwRadius(ElementType e);
  /**
   * @brief Getter for the nuclear charge.
   * @param e The ElementType.
   * @return int Returns the nuclear charge.
   */
  static int Z(ElementType e);
  /**
   * @brief Getter for the number of valence electrons.
   * @param e The ElementType.
   * @return int Returns the number of valence electrons.
   */
  static int valElectrons(ElementType e);
  /**
   * @brief Getter for the number of valence s-electrons.
   * @param e The ElementType.
   * @return int Returns the number of valence s-electrons.
   */
  static int sElectrons(ElementType e);
  /**
   * @brief Getter for the number of valence p-electrons.
   * @param e The ElementType.
   * @return int Returns the number of valence p-electrons.
   */
  static int pElectrons(ElementType e);
  /**
   * @brief Getter for the number of valence d-electrons.
   * @param e The ElementType.
   * @return int Returns the number of valence d-electrons.
   */
  static int dElectrons(ElementType e);

 private:
  /**
   * @brief A map mapping between string an enum  type of an element.
   *
   * Note: please do not use this map directly, instead use:
   * ElementInfo::elementTypeForSymbol() beacaus it includes error handling.
   */
  static std::map<std::string, ElementType> stringToElementType;
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_ELEMENTINFO_H
