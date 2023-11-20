/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ELEMENTINFO_H
#define UTILS_ELEMENTINFO_H

#include "Utils/Geometry/ElementTypes.h"
#include "Utils/MSVCCompatibility.h"
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

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
   *
   * Permissive regarding digits specifying isotope atomic mass numbers, either
   * pre- or postfixed.
   *
   * @code{cpp}
   * ElementType hydrogen = elementTypeForSymbol("H"); // ElementType::H
   * ElementType h1prefixed = elementTypeForSymbol("1H"); // ElementType::H1
   * ElementType h1postfixed = elementTypeForSymbol("H1"); // ElementType::H1
   * ElementType deuterium = elementTypeForSymbol("D"); // ElementType::D
   * ElementType h2prefixed = elementTypeForSymbol("2H"); // ElementType::D
   * ElementType tritium = elementTypeForSymbol("T"); // ElementType::T
   * ElementType h3prefixed = elementTypeForSymbol("3H"); // ElementType::T
   * @endcode
   *
   * @param symbol The atom symbol. First character captitalized, other lower case.
   * @return ElementType Returns the ElementType.
   */
  static ElementType elementTypeForSymbol(const std::string& symbol);

  /**
   * @brief Vector of all ElementType enum.
   *
   * @return std::vector<ElementType> Returns all possible ElementTypes (excluding None).
   */
  static std::vector<ElementType> allImplementedElements();

  /**
   * @brief Translation from ElementType enum to std::string.
   * @param e The ElementType.
   * @return std::string Returns the symbol string.
   */
  static std::string symbol(ElementType e);

  /*! @brief Standard atomic weight of the atom
   *
   * The standard atomic weight of an element (e.g. H) is the average of its
   * isotopic weights weighted by their natural abundance. If no natural
   * abundance for an element was measured or no standard atomic weight is
   * defined, returns the weight of one of its isotopes.
   *
   * The atomic weight of an isotope (e.g. D) is the mass of the isotope scaled
   * onto the standard atomic weight scale, where the standard atomic weight of
   * C-12 is set to 12.
   *
   * @param element The element type for which to fetch the standard atomic weight
   * @returns standard atomic weight in unified atomic mass units (u)
   */
  static double mass(ElementType e);
  /**
   * @brief Getter for the covalent radius in atomic units.
   *
   * References:
   * - "Atomic Radii of the Elements" in CRC Handbook of Chemistry and Physics, 100th Edition
   * (Internet Version 2019), John R. Rumble, ed., CRC Press/Taylor & Francis, Boca Raton, FL.
   * - DOI: 10.1039/b801115j
   * - DOI: 10.1002/chem.200800987
   * @return double  Returns the covalent radius in atomic units.
   */

  static double covalentRadius(ElementType e);
  /**
   * @brief Getter for the van der Waals radius in atomic units.
   * @param e The ElementType.
   * @return double Returns the van der Waals radius in atomic units.
   */
  static double vdwRadius(ElementType e);

  /**
   * @brief Getter for the Pauling electronegativity.
   * @param e The ElementType.
   *
   * CRC Handbook of Physics and Chemistry 85th Edition,
   * CRC Press/Taylor and Francis, Boca Raton, FL, Molecular Structure and Spectroscopy,
   * Sec. 9, page 74.\n
   * Original references:
   * Pauling, L., The Nature of the Chemical Bond, Third Edition, Cornell University Press, Ithaca, NY, 1960.\n
   * Allen, L. C., J. Am. Chem. Soc., 111, 9003, 1989.\n
   * Allred, A. L., J. Inorg. Nucl. Chem., 17, 215, 1961.
   *
   * @return double Returns the Pauling electronegativity.
   */
  static double paulingElectronegativity(ElementType e);

  /**
   * @brief Getter for the atomic number
   * @param e The ElementType.
   * @return int Returns the atomic number
   */
  static constexpr unsigned Z(const ElementType e) noexcept {
    constexpr unsigned bitmask = 127;
    return static_cast<unsigned>(e) bitand bitmask;
  }

  //! Gets the atomic mass number of an element type
  static constexpr unsigned A(const ElementType e) noexcept {
    constexpr unsigned zBits = 7;
    return static_cast<unsigned>(e) >> zBits;
  }

  /*! @brief Natural abundance of an isotope
   *
   * Note that the stored natural abundances of particular isotopes may not
   * sum to one, but may all be zero for cases in which no natural abundances
   * have been measured.
   *
   * @param element The element type of the isotope
   * @throws std::logic_error For unspecified isotopes (e.g. H, but not Be,
   *   which is monoisotopic)
   *
   * @returns natural relative abundance of an isotope in fractions <= 1
   */
  static double abundance(ElementType element);

  //! Returns the element for a particular atomic number
  static ElementType element(unsigned z);

  //! Returns an element type corresponding to atomic- and atomic mass-numbers
  static ElementType isotope(unsigned z, unsigned a);

  //! Returns the isotopes of an element, unordered
  static std::vector<ElementType> isotopes(ElementType element);

  //! Returns the base of an isotope (e.g. Li for Li6)
  static ElementType base(ElementType isotope);

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
   *
   * @note Valence electron shell assignment is done by application of
   * Madelung's rule, without incorporation of experimental exceptions.
   */
  static int sElectrons(ElementType e);

  /**
   * @brief Getter for the number of valence p-electrons.
   * @param e The ElementType.
   * @return int Returns the number of valence p-electrons.
   *
   * @note Valence electron shell assignment is done by application of
   * Madelung's rule, without incorporation of experimental exceptions.
   */
  static int pElectrons(ElementType e);

  /**
   * @brief Getter for the number of valence d-electrons.
   * @param e The ElementType.
   * @return int Returns the number of valence d-electrons.
   *
   * @note Valence electron shell assignment is done by application of
   * Madelung's rule, without incorporation of experimental exceptions.
   */
  static int dElectrons(ElementType e);

  /**
   * @brief Getter for the number of valence f-electrons.
   * @param e The ElementType.
   * @return int Returns the number of valence f-electrons.
   *
   * @note Valence electron shell assignment is done by application of
   * Madelung's rule, without incorporation of experimental exceptions.
   */
  static int fElectrons(ElementType e);

  /**
   * @brief A map mapping between string an enum  type of an element.
   *
   * @note Please do not use this map directly, instead use:
   * ElementInfo::elementTypeForSymbol() because it includes error handling.
   */
  static const std::unordered_map<std::string, ElementType>& stringToElementType();

 private:
  struct IsotopeData {
    double mass;
    double abundance;
  };

  static const std::unordered_map<ElementType, IsotopeData>& isotopeMap();
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_ELEMENTINFO_H
