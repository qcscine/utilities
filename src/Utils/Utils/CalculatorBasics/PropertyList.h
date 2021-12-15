/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_PROPERTYLIST_H
#define UTILS_PROPERTYLIST_H

#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/DataStructures/AtomicGtos.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/DipoleMatrix.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/Math/AtomicSecondDerivativeCollection.h>
#include <Utils/Math/IterativeDiagonalizer/SpinAdaptedEigenContainer.h>
#include <Utils/Properties/Thermochemistry/ThermochemistryCalculator.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Typenames.h>
#include <array>
#include <type_traits>
#include <unordered_map>

namespace Scine {
namespace Utils {

/*! @brief The properties contained are assigned a bit. This can be switched on or off
 *         to flag the presence or absence of the property. */
enum class Property : unsigned {
  Energy = (1 << 0),
  Gradients = (1 << 1),
  Hessian = (1 << 2),
  AtomicHessians = (1 << 3),
  Dipole = (1 << 4),
  DipoleGradient = (1 << 5),
  DipoleMatrixAO = (1 << 6),
  DipoleMatrixMO = (1 << 7),
  DensityMatrix = (1 << 8),
  OneElectronMatrix = (1 << 9),
  TwoElectronMatrix = (1 << 10),
  OverlapMatrix = (1 << 11),
  CoefficientMatrix = (1 << 12),
  OrbitalEnergies = (1 << 13),
  ElectronicOccupation = (1 << 14),
  Thermochemistry = (1 << 15),
  ExcitedStates = (1 << 16),
  AOtoAtomMapping = (1 << 17),
  AtomicCharges = (1 << 18),
  BondOrderMatrix = (1 << 19),
  Description = (1 << 20),
  SuccessfulCalculation = (1 << 21),
  ProgramName = (1 << 22),
  PointChargesGradients = (1 << 23),
  AtomicGtos = (1 << 24),
  GridOccupation = (1 << 25)
};

// clang-format off
using PropertyTypeTuple =
    std::tuple<
    double, /*Property::Energy*/
    GradientCollection, /*Property::Gradients*/
    HessianMatrix, /*Property::Hessian*/
    AtomicSecondDerivativeCollection, /*Property::AtomicHessians*/
    Dipole, /*Property::Dipole*/
    DipoleGradient, /*Property::DipoleGradient*/
    DipoleMatrix, /*Property::DipoleMatrixAO*/
    DipoleMatrix, /*Property::DipoleMatrixMO*/
    DensityMatrix, /*Property::DensityMatrix*/
    Eigen::MatrixXd, /*Property::OneElectronMatrix*/
    SpinAdaptedMatrix, /*Property::TwoElectronMatrix*/
    Eigen::MatrixXd, /*Property::OverlapMatrix*/
    MolecularOrbitals, /*CProperty::CoefficientMatrix*/
    SingleParticleEnergies, /*Property::OrbitalEnergies*/
    LcaoUtils::ElectronicOccupation, /*Property::ElectronicOccupation*/
    ThermochemicalComponentsContainer, /*Property::Thermochemistry*/
    SpinAdaptedElectronicTransitionResult, /*Property::ExcitedStates*/
    AtomsOrbitalsIndexes, /*Property::AOtoAtomMapping*/
    std::vector<double>, /*Property::AtomicCharges*/
    BondOrderCollection, /*Property::BondOrderMatrix*/
    std::string, /*Property::Description*/
    bool, /*Property::SuccessfulCalculation*/
    std::string, /*Property::ProgramName*/
    GradientCollection, /*Property::PointChargesGradients*/
    std::unordered_map<int, AtomicGtos>, /*Property::AtomicGtos*/
    std::vector<int> /*Property::GridOccupation*/
    >;
// clang-format on

static_assert(std::tuple_size<PropertyTypeTuple>::value == 26,
              "Tuple does not contain as many elements as there are properties");

constexpr std::array<Property, std::tuple_size<PropertyTypeTuple>::value> allProperties{{Property::Energy,
                                                                                         Property::Gradients,
                                                                                         Property::Hessian,
                                                                                         Property::AtomicHessians,
                                                                                         Property::Dipole,
                                                                                         Property::DipoleGradient,
                                                                                         Property::DipoleMatrixAO,
                                                                                         Property::DipoleMatrixMO,
                                                                                         Property::DensityMatrix,
                                                                                         Property::OneElectronMatrix,
                                                                                         Property::TwoElectronMatrix,
                                                                                         Property::OverlapMatrix,
                                                                                         Property::CoefficientMatrix,
                                                                                         Property::OrbitalEnergies,
                                                                                         Property::ElectronicOccupation,
                                                                                         Property::Thermochemistry,
                                                                                         Property::ExcitedStates,
                                                                                         Property::AOtoAtomMapping,
                                                                                         Property::AtomicCharges,
                                                                                         Property::BondOrderMatrix,
                                                                                         Property::Description,
                                                                                         Property::SuccessfulCalculation,
                                                                                         Property::ProgramName,
                                                                                         Property::PointChargesGradients,
                                                                                         Property::AtomicGtos,
                                                                                         Property::GridOccupation}};

// Python binding names
constexpr std::array<const char*, std::tuple_size<PropertyTypeTuple>::value> allPropertyNames{"energy",
                                                                                              "gradients",
                                                                                              "hessian",
                                                                                              "atomic_hessian",
                                                                                              "dipole",
                                                                                              "dipole_gradient",
                                                                                              "ao_dipole_matrix",
                                                                                              "mo_dipole_matrix",
                                                                                              "density_matrix",
                                                                                              "one_electron_matrix",
                                                                                              "two_electron_matrix",
                                                                                              "overlap_matrix",
                                                                                              "coefficient_matrix",
                                                                                              "orbital_energies",
                                                                                              "electronic_occupation",
                                                                                              "thermochemistry",
                                                                                              "excited_states",
                                                                                              "ao_to_atom_mapping",
                                                                                              "atomic_charges",
                                                                                              "bond_orders",
                                                                                              "description",
                                                                                              "successful_calculation",
                                                                                              "program_name",
                                                                                              "point_charges_gradients",
                                                                                              "atomic_gtos",
                                                                                              "grid_occupation"};

/* other variants of doing this:
 * - Use a constexpr map datatype
 * - Use a std::array<std::pair<Property, const char*>> and linear-search it
 * - Use a separate array of strings and figure out the index from the original property
 */

constexpr unsigned getPropertyIndex(Property property) {
  unsigned index = allProperties.size();
  for (unsigned i = 0; i < allProperties.size(); ++i) {
    if (allProperties.at(i) == property) {
      index = i;
      break;
    }
  }

  if (index == allProperties.size()) {
    throw std::logic_error("constexpr failed to find property" +
                           std::to_string(static_cast<std::underlying_type<Property>::type>(property)));
  }

  return index;
}

static_assert(getPropertyIndex(Property::Gradients) == 1, "Fn doesn't work");

constexpr const char* propertyTypeName(Property property) {
  unsigned enumIndex = getPropertyIndex(property);
  return allPropertyNames.at(enumIndex);
}

template<Property property>
struct PropertyType {
  using type = std::tuple_element_t<getPropertyIndex(property), PropertyTypeTuple>;
  using Type = type;
  static constexpr const char* name = propertyTypeName(property);
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
