/**
 * @file ElementData.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ELEMENTDATA_H_
#define UTILS_ELEMENTDATA_H_

#include "Utils/Constants.h"
#include "Utils/Geometry/ElementTypes.h"
#include <map>
#include <memory>
#include <string>
#include <utility>

namespace Scine {
namespace Utils {
namespace Constants {

/**
 * @class ElementDataSingleton ElementData.h
 * @brief Provides a mapping of ElementType to ElementData and accessing Elements by a symbol string.
 *
 * This class is a singleton to avoid multiple instances of the same hardcoded data in memory.
 *
 * Features:
 *  - Fast lookup for subscript operator [ElementType]. Throws std::out_of_range.
 *  - Slow lookup for subscript operator [std::string]. Throws ElementData::DataNotAvailable.
 */
class ElementDataSingleton {
 public:
  /**
   * @class ElementData ElementData.h
   * @brief Data type for each element of the periodic table, Data includes: Mass (a.u.), Z,
   */
  class ElementData {
   public:
    /**
     * @class DataNotAvailable ElementData.h
     * @brief An exception if the data does not exist.
     */
    class DataNotAvailable {};

    /**
     * @brief Constructor for default (empty) element.
     */
    ElementData() = default;

    /**
     * @brief Construct a new ElementData object.
     *
     * @param symbol The element symbol.
     * @param Z      The nuclear charge.
     * @param mass   The mass (isotope average, precision: 3 digits).
     * @param covalentRadiusInPicometers The covalent radius in pm.
     * @param vdWRadiusInPicometers The van der Waals radius in pm.
     * @param valElectrons The number of electrons in the valence shell.
     * @param sElectrons The number of s-electrons in the valence shell.
     * @param pElectrons The number of p-electrons in the valence shell.
     * @param dElectrons The number of d-electrons in the valence shell.
     */
    ElementData(std::string symbol, int Z, double mass, double covalentRadiusInPicometers = -1,
                double vdWRadiusInPicometers = -1, int valElectrons = -1, int sElectrons = -1, int pElectrons = -1,
                int dElectrons = -1)
      : d_symbol(std::move(symbol)),
        d_Z(Z),
        d_mass(mass),
        d_covalentRadius(covalentRadiusInPicometers / 100 * bohr_per_angstrom),
        d_vdWRadius(vdWRadiusInPicometers / 100 * bohr_per_angstrom),
        d_valElectrons(valElectrons),
        d_sElectrons(sElectrons),
        d_pElectrons(pElectrons),
        d_dElectrons(dElectrons){};

    /**
     * @brief Getter for the element symbol.
     * @return std::string Returns the element symbol as string.
     */
    std::string symbol() const {
      return d_symbol;
    }
    /**
     * @brief Getter for the atomic number.
     * @return int Returns the atomic number.
     */
    int Z() const {
      return d_Z;
    }
    /**
     * @brief Getter for the mass in atomic mass units (u).
     * @return double Returns the atomic mass (isotope average, precision: 3 digits)
     */
    double mass() const {
      return d_mass;
    }
    /**
     * @brief Getter for the covalent radius.
     *
     * References:
     * - "Atomic Radii of the Elements" in CRC Handbook of Chemistry and Physics, 100th Edition
     * (Internet Version 2019), John R. Rumble, ed., CRC Press/Taylor & Francis, Boca Raton, FL.
     * - DOI: 10.1039/b801115j
     * - DOI: 10.1002/chem.200800987
     * @return double  Returns the covalent radius in atomic units.
     */
    double covalentRadius() const {
      if (d_covalentRadius > 0.0)
        return d_covalentRadius;
      throw DataNotAvailable();
    }
    /**
     * @brief Getter for the Van-der-Waals Radius.
     * @return double  Returns the Van-der-Waals Radius in atomic units.
     */
    double vdWRadius() const {
      if (d_vdWRadius > 0.0)
        return d_vdWRadius;
      throw DataNotAvailable();
    }
    /**
     * @brief Getter for the number of valence electrons.
     * @return int  Returns the number of valence electrons.
     */
    int valElectrons() const {
      if (d_valElectrons > -1)
        return d_valElectrons;
      throw DataNotAvailable();
    }
    /**
     * @brief Getter for the number of valence s-electrons.
     * @return int  Returns the number of valence s-electrons.
     */
    int sElectrons() const {
      if (d_sElectrons > -1)
        return d_sElectrons;
      throw DataNotAvailable();
    }
    /**
     * @brief Getter for the number of valence p-electrons.
     * @return int  Returns the number of valence p-electrons.
     */
    int pElectrons() const {
      if (d_pElectrons > -1)
        return d_pElectrons;
      throw DataNotAvailable();
    }
    /**
     * @brief Getter for the number of valence d-electrons.
     * @return int  Returns the number of valence d-electrons.
     */
    int dElectrons() const {
      if (d_dElectrons > -1)
        return d_dElectrons;
      throw DataNotAvailable();
    }

   private:
    std::string d_symbol{};

    int d_Z{0};
    double d_mass{-1};
    double d_covalentRadius{-1};
    double d_vdWRadius{-1};

    int d_valElectrons{-1};
    int d_sElectrons{-1};
    int d_pElectrons{-1};
    int d_dElectrons{-1};
  };

 public:
  /**
   * @brief Access singleton instance. Creates one if necessary.
   * @return const ElementDataSingleton& A reference to the one instance allowed.
   */
  static const ElementDataSingleton& instance();
  /**
   * @brief Access element information based on element symbol.
   * @param symbol The symbol of the element (1st character upper case, 2nd lower case).
   * @throws ElementSymbolNotFound

   * @return const ElementData& Returns the data associated with this element.
   */
  const ElementData& operator[](const std::string& symbol) const;
  /**
   * @brief Access element information based on type. Fastest lookup.
   * @param type The ElementType.
   * @throws std::out_of_range
   * @return const ElementData& Returns the data associated with this element.
   */
  const ElementData& operator[](const ElementType& type) const;

 private:
  /// @brief Private constructor to prevent instantiation by client.
  ElementDataSingleton();

  /// @brief The singleton instance.
  static std::unique_ptr<ElementDataSingleton> d_instance;

  /// @brief Storage: internal map.
  std::map<ElementType, ElementData> d_container;

  /// @brief Creates hard coded element data.
  void init_data();
};

} /* namespace Constants */
} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_ELEMENTDATA_H_
