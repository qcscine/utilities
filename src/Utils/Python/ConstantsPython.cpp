/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/ElementData.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_constants(pybind11::module& m) {
  m.attr("ELEMENTARY_CHARGE") = Constants::elementaryCharge;
  m.attr("AVOGADRO_NUMBER") = Constants::avogadroNumber;
  m.attr("PI") = Constants::pi;
  m.attr("ATOMIC_MASS_UNIT") = Constants::atomicMassUnit;
  m.attr("ELECTRON_REST_MASS") = Constants::electronRestMass;
  m.attr("INVERSE_FINE_STRUCTURE_CONSTANT") = Constants::inverseFineStructureConstant;
  m.attr("RAD_PER_DEGREE") = Constants::rad_per_degree;
  m.attr("DEGREE_PER_RAD") = Constants::degree_per_rad;
  m.attr("METER_PER_BOHR") = Constants::meter_per_bohr;
  m.attr("BOHR_PER_METER") = Constants::bohr_per_meter;
  m.attr("ANGSTROM_PER_METER") = Constants::angstrom_per_meter;
  m.attr("METER_PER_ANGSTROM") = Constants::meter_per_angstrom;
  m.attr("ANGSTROM_PER_BOHR") = Constants::angstrom_per_bohr;
  m.attr("BOHR_PER_ANGSTROM") = Constants::bohr_per_angstrom;
  m.attr("HARTREE_PER_EV") = Constants::hartree_per_ev;
  m.attr("EV_PER_HARTREE") = Constants::ev_per_hartree;
  m.attr("JOULE_PER_HARTREE") = Constants::joule_per_hartree;
  m.attr("HARTREE_PER_JOULE") = Constants::hartree_per_joule;
  m.attr("JOULE_PER_CALORIE") = Constants::joule_per_calorie;
  m.attr("CALORIE_PER_JOULE") = Constants::calorie_per_joule;
  m.attr("INVERSE_CENTIMETER_PER_HARTREE") = Constants::invCentimeter_per_hartree;
  m.attr("HARTREE_PER_INVERSE_CENTIMETER") = Constants::hartree_per_invCentimeter;
  m.attr("KJPERMOL_PER_HARTREE") = Constants::kJPerMol_per_hartree;
  m.attr("HARTREE_PER_KJPERMOL") = Constants::hartree_per_kJPerMol;
  m.attr("KCALPERMOL_PER_HARTREE") = Constants::kCalPerMol_per_hartree;
  m.attr("HARTREE_PER_KCALPERMOL") = Constants::hartree_per_kCalPerMol;
  m.attr("U_PER_KG") = Constants::u_per_kg;
  m.attr("KG_PER_U") = Constants::kg_per_u;
  m.attr("ELECTRONRESTMASS_PER_KG") = Constants::electronRestMass_per_kg;
  m.attr("KG_PER_ELECTRONRESTMASS") = Constants::kg_per_electronRestMass;
  m.attr("ELECTRONRESTMASS_PER_U") = Constants::electronRestMass_per_u;
  m.attr("U_PER_ELECTRONRESTMASS") = Constants::u_per_electronRestMass;
}
