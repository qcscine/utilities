/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 * @brief This header file defines constants commonly used in computational chemistry.
 *
 * Source: http://physics.nist.gov/cuu/Constants/Table/allascii.txt, 03.11.2015.
 */
#ifndef UTILS_CONSTANTS_H_
#define UTILS_CONSTANTS_H_

#include "Utils/Technical/StrongType.h"

namespace Scine {
namespace Utils {
/**
 * @brief A namespace for all constant (hardcoded) data
 *
 * This namespace does not include fitted parameters for
 * specific methods, only general constant parameters,
 * such as natural constants and atomic data.
 */
namespace Constants {
/*=====================
 *  Natural Constants
 *=====================*/
// Constants
constexpr double elementaryCharge = 1.6021766208e-19; // C
constexpr double avogadroNumber = 6.022140857e23;     // mol^-1
constexpr double bohrRadius = 0.52917721067e-10;      // m
constexpr double pi = 3.14159265358979323846;

/*======================
 *  Conversion Factors
 *======================*/
// Rad and Degree
constexpr double rad_per_degree = pi / 180;
constexpr double degree_per_rad = 180 / pi;

// Length
constexpr double meter_per_bohr = bohrRadius;
constexpr double bohr_per_meter = 1 / meter_per_bohr;

constexpr double angstrom_per_meter = 1e10;
constexpr double meter_per_angstrom = 1 / angstrom_per_meter;

constexpr double angstrom_per_bohr = angstrom_per_meter * meter_per_bohr;
constexpr double bohr_per_angstrom = 1. / angstrom_per_bohr;

// Energy
constexpr double hartree_per_ev = 3.674932248e-2;
constexpr double ev_per_hartree = 1 / hartree_per_ev;

constexpr double joule_per_hartree = 4.359744650e-18;
constexpr double hartree_per_joule = 1 / joule_per_hartree;

constexpr double joule_per_calorie = 4.184;
constexpr double calorie_per_joule = 1 / joule_per_calorie;

constexpr double kJPerMol_per_hartree = joule_per_hartree / 1000 * avogadroNumber;
constexpr double hartree_per_kJPerMol = 1 / kJPerMol_per_hartree;

constexpr double kCalPerMol_per_hartree = joule_per_hartree * calorie_per_joule / 1000 * avogadroNumber;
constexpr double hartree_per_kCalPerMol = 1 / kCalPerMol_per_hartree;

// Spectroscopic conversions E = hcw
constexpr double ev_per_invCentimeter = 8.065540107e3;
constexpr double invCentimeter_per_ev = 1.239842573e-4;

constexpr double hartree_per_invCentimeter = 4.556335281e-6;
constexpr double invCentimeter_per_hartree = 219474.63;

constexpr double kJPerMol_per_invCentimeter = 2.859144166e-3;
constexpr double invCentimeter_per_kJPerMol = 349.7550112241469;

constexpr double kCalPerMol_per_invCentimeter = 1.196265919e-2;
constexpr double invCentimeter_per_kCalPerMol = 8.359345393;

} /* namespace Constants */
/*=========
 *  Units
 *=========*/
// Length
using Bohr = StrongType<double, struct BohrType>;
using Angstrom = StrongType<double, struct AngstromType>;

// Energy
using Hartree = StrongType<double, struct HartreeType>;
using kJPerMol = StrongType<double, struct kJPerMolType>;
using kCalPerMol = StrongType<double, struct kCalPerMolType>;
using eV = StrongType<double, struct eVType>;

/*========================
 *  Conversion Functions
 *========================*/

// Length
/**
 * @brief Conversion from Bohr to Angstrom.
 * @deprecated
 * @param b Length in/as Bohr.
 * @return constexpr double Returns the length in Angstrom.
 */
constexpr double toAngstrom(Bohr b) {
  return b.get() * Constants::angstrom_per_bohr;
}
/**
 * @brief Conversion from Angstrom to Bohr.
 * @deprecated
 * @param a Length in/as Angstrom.
 * @return constexpr double Returns the length in Bohr.
 */
constexpr double toBohr(Angstrom a) {
  return a.get() * Constants::bohr_per_angstrom;
}

// Energy
/**
 * @brief Conversion from Hartree to kJ/mol.
 * @deprecated
 * @param h Energy in/as Hartree.
 * @return constexpr double Returns the energy in kJ/mol.
 */
constexpr double toKJPerMol(Hartree h) {
  return h.get() * Constants::kJPerMol_per_hartree;
}
/**
 * @brief Conversion from kJ/mol to Hartree.
 * @deprecated
 * @param kjpm Energy in/as kJ/mol.
 * @return constexpr double Returns the energy in Hartree.
 */
constexpr double toHartree(kJPerMol kjpm) {
  return kjpm.get() * Constants::hartree_per_kJPerMol;
}
/**
 * @brief Conversion from Hartree to kcal/mol.
 * @deprecated
 * @param h Energy in/as Hartree.
 * @return constexpr double Returns the energy in kcal/mol.
 */
constexpr double toKCalPerMol(Hartree h) {
  return h.get() * Constants::kCalPerMol_per_hartree;
}
/**
 * @brief Conversion from kcal/mol to Hartree.
 * @deprecated
 * @param kcpm Energy in/as kcal/mol.
 * @return constexpr double Returns the energy in Hartree.
 */
constexpr double toHartree(kCalPerMol kcpm) {
  return kcpm.get() * Constants::hartree_per_kCalPerMol;
}

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_CONSTANTS_H_
