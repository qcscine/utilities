/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_PROPERTYLIST_H
#define UTILS_PROPERTYLIST_H

#include <type_traits>

namespace Scine {
namespace Utils {

/*! @brief The properties contained are assigned a bit. This can be switched on or off
 *         to flag the presence or absence of the property. */
enum class Property {
  Energy = (1 << 0),
  Gradients = (1 << 1),
  Hessian = (1 << 2),
  Dipole = (1 << 3),
  DipoleGradient = (1 << 4),
  DipoleMatrixAO = (1 << 5),
  DipoleMatrixMO = (1 << 6),
  OneElectronMatrix = (1 << 7),
  TwoElectronMatrix = (1 << 8),
  BondOrderMatrix = (1 << 9)
};

/*! @brief Returns a Property object that is the superset of the two properties given as argument*/
constexpr inline Property operator|(Property v1, Property v2);
/*! @brief Returns a Property object that is the subset of the two properties given as argument*/
constexpr inline bool operator&(Property v1, Property v2);

/*!
 * This class defines a list of properties that can be calculated in a single-point calculation.
 */
class PropertyList {
 public:
  //! Initializes the property enum to be empty, i.e. all switched off
  PropertyList() : properties_(static_cast<Property>(0)) {
  }

  /*! Constructor from properties; not explicit to allow for automatic conversion. */
  PropertyList(Property p) : properties_(p) {
  }

  /*! Checks for the presence of a PropertyList given as argument as subset of the current object. */
  bool containsSubSet(const PropertyList& pl) const {
    auto combinedProperties = pl.properties_ | properties_;
    return combinedProperties == properties_;
  }

  /*! Switches on the bits that are switched on in the argument Property v */
  void addProperty(const Property v) {
    properties_ = properties_ | v;
  }

 private:
  Property properties_;
};

/*! Allow combination of properties. */
constexpr inline Property operator|(const Property v1, const Property v2) {
  using utype = std::underlying_type<Property>::type;
  return static_cast<Property>(static_cast<utype>(v1) | static_cast<utype>(v2));
}

/*! Allow to check if there is a flag overlap. */
constexpr inline bool operator&(const Property v1, const Property v2) {
  using utype = std::underlying_type<Property>::type;
  return (static_cast<utype>(v1) & static_cast<utype>(v2)) != 0;
}

} // namespace Utils
} // namespace Scine
#endif // UTILS_PROPERTYLIST_H
