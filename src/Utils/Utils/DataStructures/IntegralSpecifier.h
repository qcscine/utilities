/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_INTEGRALSPECIFIER_H
#define UTILSOS_INTEGRALSPECIFIER_H

/* external */
#include <boost/optional.hpp>
#include <unordered_map>
/* internal */
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"

namespace Scine {
namespace Utils {
namespace Integrals {

/**
 * @brief `first` of ReturnKey. Is used to specify the matrix returned by the integral routine.
 *         `none` is for scalar integrals. `x, y, z` for vector-integrals, and so on.
 */
enum class Component { none = 0, x = 1, y = 2, z = 3, xx = 4, xy = 5, xz = 6, yy = 7, yz = 8, zz = 9 };

/**
 * @brief `second` of ReturnKey. Is used to specify the matrix returned by the integral routine.
 *        `value` if no derivative order is 0, `x1, ...` for derivative order 1.
 */
enum class DerivKey { value = 0, x = 0, y = 1, z = 2 };

using Center = std::size_t;

using ReturnKey = std::tuple<Component, DerivKey, Center>;

using Spin = int;

/**
 * Particle type struct
 * Only relevant for pre-Born--Oppenheimer program.
 * In standard cases, the particle type is electron by default.
 *
 * For input files:
 * Always provide the `symbol`, which is defines in capital letter followed by small letter of element type,
 * with N as the isotope number.
 * XyN
 * Electron is E
 * Positron is Ps
 */
struct ParticleType {
  int charge;
  Spin spin; // Spin is always given as spin*2
  ElementType symbol;
  double mass;
};

/**
 * @brief returns a vector of all available `ParticleType` objects.
 * @return vector of ParticleType objects
 */
static const std::vector<ParticleType>& getParticleTypeInfo() {
  static std::vector<ParticleType> particleTypeInfo{
      // The electron is an element.
      {-1, 1, ElementType::E, 1},
      // {1, 1, Symbol::Ps, 1},
      {1, 1, ElementType::H, 1836.152701},        //  from Eur. Phys. J. D, 12, 449-466 (2000)
      {1, 2, ElementType::D, 3670.483014},        //  CODATA 2014
      {1, 1, ElementType::T, 5496.9215358826},    //  CODATA 2014
      {2, 0, ElementType::He4, 7294.29951423954}, //  CODATA 2005 - Eq. (9)
      // Nucleus mass in atomic units: atomic masss / factor - number of electrons * electron mass, i.e, charge
      {3, 3, ElementType::Li, 7.016003436645 / 0.000548579911112 - 3},    // PHYSICS.NIST.GOV
      {4, 3, ElementType::Be9, 9.01218306582 / 0.000548579911112 - 4},    // PHYSICS.NIST.GOV
      {5, 3, ElementType::B11, 11.00930516713 / 0.000548579911112 - 5},   // PHYSICS.NIST.GOV
      {6, 0, ElementType::C12, 21868.663768752986},                       // CODATA 2005
      {6, 1, ElementType::C13, 13.0033548350723 / 0.000548579911112 - 6}, // PHYSICS.NIST.GOV
      {7, 1, ElementType::N15, 15.0001088988864 / 0.000548579911112 - 7}, // PHYSICS.NIST.GOV
      {8, 1, ElementType::O, 15.0030656253 / 0.000548579911112 - 8},      // PHYSICS.NIST.GOV
      {8, 0, ElementType::O16, 29148.949593448239},                       // CODATA 2005
      {9, 1, ElementType::F19, 18.9984031627392 / 0.000548579911112 - 9}  // PHYSICS.NIST.GOV
  };
  return particleTypeInfo;
}

/**
 * Converts a string to a `ParticleType` `Symbol`
 * @param str_symbol
 * @return `Symbol`
 */
inline auto getParticleType(std::string str_symbol) -> ParticleType {
  // Make all characters of the symbol lowercase if it isn't yet
  std::transform(std::begin(str_symbol), std::end(str_symbol), std::begin(str_symbol),
                 [](const auto c) { return std::tolower(c); });
  auto str2symb = ElementInfo::stringToElementType();
  auto symb = str2symb[str_symbol];

  auto particleTypeInfo = getParticleTypeInfo();

  for (const auto& tp : particleTypeInfo) {
    if (tp.symbol == symb) {
      return tp;
    }
  }
  throw std::runtime_error("Particle type " + ElementInfo::symbol(symb) + " not implemented.");
}

/*
 * Standard integrals:
 * Overlap, Kinetic, PointCharges, Coulomb, Dipole
 * Pre-BO specific:
 * KineticCOM   Subtracts the total momentum for the elimination of the COM motion.
 * CoulombCOM                 -"-
 *              Note: those integrals do not have 8-fold symmetry!
 */
enum class Operator { Overlap, Kinetic, PointCharges, Coulomb, KineticCOM, CoulombCOM, Dipole };

/**
 * @struct IntegralSpecifier IntegralSpecifier.h
 * @brief Is passed to the Integral Interface and contains all relevant infos to evaluate the required integral.
 */
struct IntegralSpecifier {
  // Operator specifies the type integral.
  Operator op;
  // By default: use electron as paricle type.
  std::vector<ParticleType> typeVector{{-1, 1, ElementType::E, 1}};
  // Derivative order. Supported are only up to first order derivatives.
  std::size_t derivOrder = 0;
  // Optional: molecular geometry for point charges integrals.
  boost::optional<AtomCollection> atoms;
  // Optional: total mass for pre-Born--Oppenheimer integrals only.
  boost::optional<double> totalMass;
  // Optional: origin for multipole integrals.
  boost::optional<Position> multipoleOrigin;
};

} // namespace Integrals
} // namespace Utils
} // namespace Scine

#endif // UTILSOS_INTEGRALSPECIFIER_H
