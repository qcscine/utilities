/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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
#include <Utils/DataStructures/PartialHessian.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/ExternalQC/Orca/MoessbauerParameters.h>
#include <Utils/Math/AtomicSecondDerivativeCollection.h>
#include <Utils/Math/IterativeDiagonalizer/SpinAdaptedEigenContainer.h>
#include <Utils/Properties/Thermochemistry/ThermochemistryCalculator.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <array>
#include <type_traits>
#include <unordered_map>

namespace Scine {
namespace Utils {

/*! @brief The properties contained are assigned a bit. This can be switched on or off
 *         to flag the presence or absence of the property.
 *
 *  The total number of possible properties can be increased from 64 by changing the data
 *  type to a wrapped implementation of a std::bitset, see https://fschoenberger.dev/enum-bitset/
 *  for an example. Note that you have to be careful with the bit shift operations here.
 *  If you do (1 << 31) you get an integer overflow because 1 is only a unsigned int. You have to
 *  do (1ULL << 31) to ensure that 1 is actually an unsigned long long.
 *
 *  See also ResultsPython.cpp for more bit shift operations.
 *
 */
enum class Property : unsigned long long {
  Energy = (1ULL << 0),
  Gradients = (1ULL << 1),
  Hessian = (1ULL << 2),
  PartialHessian = (1ULL << 3),
  AtomicHessians = (1ULL << 4),
  Dipole = (1ULL << 5),
  DipoleGradient = (1ULL << 6),
  DipoleMatrixAO = (1ULL << 7),
  DipoleMatrixMO = (1ULL << 8),
  DensityMatrix = (1ULL << 9),
  OneElectronMatrix = (1ULL << 10),
  TwoElectronMatrix = (1ULL << 11),
  OverlapMatrix = (1ULL << 12),
  CoefficientMatrix = (1ULL << 13),
  OrbitalEnergies = (1ULL << 14),
  ElectronicOccupation = (1ULL << 15),
  Thermochemistry = (1ULL << 16),
  ExcitedStates = (1ULL << 17),
  AOtoAtomMapping = (1ULL << 18),
  AtomicCharges = (1ULL << 19),
  BondOrderMatrix = (1ULL << 20),
  Description = (1ULL << 21),
  SuccessfulCalculation = (1ULL << 22),
  ProgramName = (1ULL << 23),
  PointChargesGradients = (1ULL << 24),
  AtomicGtos = (1ULL << 25),
  GridOccupation = (1ULL << 26),
  StressTensor = (1ULL << 27),
  MoessbauerParameter = (1ULL << 28),
  PartialEnergies = (1ULL << 29),
  PartialGradients = (1ULL << 30),
  OrbitalFragmentPopulations = (1ULL << 31)
};

constexpr int N_PROPERTIES = 32;

// clang-format off
using PropertyTypeTuple =
    std::tuple<
    double, /*Property::Energy*/
    GradientCollection, /*Property::Gradients*/
    HessianMatrix, /*Property::Hessian*/
    PartialHessian, /*Property::PartialHessian*/
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
    std::vector<int>, /*Property::GridOccupation*/
    Eigen::Matrix3d, /*Property::StressTensor*/
    ExternalQC::Moessbauer::MoessbauerParameterContainer, /*Property::MoessbauerParameter*/
    std::unordered_map<std::string, double>, /*Property::PartialEnergies*/
    std::unordered_map<std::string, GradientCollection>, /*Property::PartialGradients*/
    SpinAdaptedMatrix /*Property::OrbitalFragmentPopulations*/
    >;
// clang-format on

static_assert(std::tuple_size<PropertyTypeTuple>::value == N_PROPERTIES,
              "Tuple does not contain as many elements as there are properties");

constexpr std::array<Property, N_PROPERTIES> allProperties{{Property::Energy,
                                                            Property::Gradients,
                                                            Property::Hessian,
                                                            Property::PartialHessian,
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
                                                            Property::GridOccupation,
                                                            Property::StressTensor,
                                                            Property::MoessbauerParameter,
                                                            Property::PartialEnergies,
                                                            Property::PartialGradients,
                                                            Property::OrbitalFragmentPopulations}};

static_assert(std::tuple_size<decltype(allProperties)>::value == N_PROPERTIES,
              "allProperties does not contain as many elements as there are properties");

// Python binding names
constexpr std::array<const char*, std::tuple_size<PropertyTypeTuple>::value> allPropertyNames{
    "energy",
    "gradients",
    "hessian",
    "partial_hessian",
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
    "grid_occupation",
    "stress_tensor",
    "moessbauer_parameter",
    "partial_energies",
    "partial_gradients",
    "orbital_fragment_populations"};

static_assert(std::tuple_size<decltype(allPropertyNames)>::value == N_PROPERTIES,
              "allPropertyNames does not contain as many elements as there are properties");

/* other variants of doing this:
 * - Use a constexpr map datatype
 * - Use a std::array<std::pair<Property, const char*>> and linear-search it
 * - Use a separate array of strings and figure out the index from the original property
 */

constexpr unsigned int getPropertyIndex(Property property) {
  unsigned index = allProperties.size();
  for (unsigned i = 0; i < allProperties.size(); ++i) {
    if (allProperties.at(i) == property) {
      index = i;
      break;
    }
  }

  if (index == allProperties.size()) {
    throw std::logic_error("constexpr failed to find property " +
                           std::to_string(static_cast<std::underlying_type<Property>::type>(property)));
  }

  return index;
}

static_assert(getPropertyIndex(Property::Gradients) == 1, "Fn doesn't work");
static_assert(getPropertyIndex(Property::OrbitalFragmentPopulations) == 31, "Fn doesn't work");

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
/*! @brief Returns a Property object that is the XOR combination of the two properties given as argument*/
constexpr inline Property operator^(Property v1, Property v2);

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

  PropertyList intersection(const PropertyList& pl) const {
    // & operator is already overloaded with a bool return type
    // so we duplicate the code here
    using utype = std::underlying_type<Property>::type;
    return static_cast<Property>(static_cast<utype>(properties_) & static_cast<utype>(pl.properties_));
  }

  /*! Switches on the bits that are switched on in the argument Property v */
  void addProperty(const Property v) {
    properties_ = properties_ | v;
  }
  /*! Switches on the bits that are switched on in the argument Property v */
  void addProperties(const PropertyList& v) {
    properties_ = properties_ | v.properties_;
  }
  /*! Switches off the bits that are switched on in the argument Property v */
  void removeProperty(const Property v) {
    if (!containsSubSet(v)) {
      // not present anyway
      return;
    }
    properties_ = properties_ ^ v;
  }
  void removeProperties(const PropertyList v) {
    if (!containsSubSet(v)) {
      // not present anyway
      return;
    }
    properties_ = properties_ ^ v.properties_;
  }

  std::string print() const {
    std::string p = "[";
    for (const auto& prop : allProperties) {
      if (this->containsSubSet(prop)) {
        if (p.size() > 1) {
          p += ", ";
        }
        p += std::string(propertyTypeName(prop));
      }
    }
    return p + "]";
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 private:
  Property properties_;
};

/*! Allow combination of properties. */
constexpr inline Property operator|(const Property v1, const Property v2) {
  using utype = std::underlying_type<Property>::type;
  return static_cast<Property>(static_cast<utype>(v1) | static_cast<utype>(v2));
}

/*! Allows to eliminate a property based on XOR logic. */
constexpr inline Property operator^(const Property v1, const Property v2) {
  using utype = std::underlying_type<Property>::type;
  return static_cast<Property>(static_cast<utype>(v1) ^ static_cast<utype>(v2));
}

/*! Allow to check if there is a flag overlap. */
constexpr inline bool operator&(const Property v1, const Property v2) {
  using utype = std::underlying_type<Property>::type;
  return (static_cast<utype>(v1) & static_cast<utype>(v2)) != 0;
}

} // namespace Utils
} // namespace Scine
#endif // UTILS_PROPERTYLIST_H
