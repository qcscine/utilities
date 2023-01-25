/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Bonds/BondDetector.h>
#include <Utils/Bonds/SolidStateBondDetector.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_bond_detector(pybind11::module& m) {
  pybind11::class_<BondDetector> bond_detector(m, "BondDetector",
                                               R"delim(
      A class to detect bonds based on interatomic distances and covalent radii.

      A bond is detected if the distance between two atoms is smaller than the sum of their covalent radii plus 0.4 Angstrom.
      A binary decision on whether a bond exists (resulting in a bond order of 1.0) or not (yielding a bond order of 0.0) is made.

      If periodic boundaries are considered, bonds across these boundaries can optionally be given the order of -1.0.

      The covalent radii were extracted from the Cambridge Structural Database (CSD) on 04/08/2020:
      https://www.ccdc.cam.ac.uk/support-and-resources/ccdcresources/Elemental_Radii.xlsx

      All method calls also allow to detect bonds based on van der Waals radii based on an optional flag.
      The van der Waals radii are specified in scine_utilities.ElementInfo

      References:
      E. C. Meng, R. A. Lewis, Comput. Chem. 1991, 12, 891-898. [DOI: 10.1002/jcc.540120716]
      C. R. Groom, I. J. Bruno, M. P. Lightfoot and S. C. Ward, Acta Cryst. 2016, B72, 171-179. [DOI: 10.1107/S2052520616003954]
    )delim");
  bond_detector.def_static("detect_bonds", pybind11::overload_cast<const AtomCollection&, bool>(&BondDetector::detectBonds),
                           pybind11::arg("atom_collection"), pybind11::arg("van_der_waals_bond") = false, R"delim(
      Generates a BondOrderCollection from an AtomCollection based on interatomic distances.)delim");
  bond_detector.def_static("detect_bonds",
                           pybind11::overload_cast<const PeriodicSystem&, bool, bool>(&BondDetector::detectBonds),
                           pybind11::arg("periodic_system"), pybind11::arg("bonds_across_boundaries_negative") = false,
                           pybind11::arg("van_der_waals_bond") = false, R"delim(
      Generates a BondOrderCollection from a PeriodicSystem based on interatomic distances and periodic boundary conditions.)delim");
  bond_detector.def_static(
      "detect_bonds",
      pybind11::overload_cast<const AtomCollection&, const PeriodicBoundaries&, bool, bool>(&BondDetector::detectBonds),
      pybind11::arg("atom_collection"), pybind11::arg("pbc"), pybind11::arg("bonds_across_boundaries_negative") = false,
      pybind11::arg("van_der_waals_bond") = false, R"delim(
      Generates a BondOrderCollection from an AtomCollection based on interatomic distances and periodic boundary conditions.)delim");
  bond_detector.def_static(
      "detect_bonds",
      pybind11::overload_cast<const ElementTypeCollection&, const PositionCollection&, bool>(&BondDetector::detectBonds),
      pybind11::arg("elements"), pybind11::arg("positions"), pybind11::arg("van_der_waals_bond") = false, R"delim(
      Generates a BondOrderCollection from an ElementTypeCollection and a PositionCollection based on interatomic distances.)delim");
  bond_detector.def_static(
      "detect_bonds",
      pybind11::overload_cast<const ElementTypeCollection&, const PositionCollection&, const PeriodicBoundaries&, bool, bool>(
          &BondDetector::detectBonds),
      pybind11::arg("elements"), pybind11::arg("positions"), pybind11::arg("pbc"),
      pybind11::arg("bonds_across_boundaries_negative") = false, pybind11::arg("van_der_waals_bond") = false, R"delim(
      Generates a BondOrderCollection from an ElementTypeCollection and a PositionCollection based on interatomic distances and periodic boundary conditions.)delim");
  bond_detector.def_static(
      "bond_exists",
      pybind11::overload_cast<const ElementType&, const ElementType&, const Position&, const Position&, bool>(
          &BondDetector::bondExists),
      pybind11::arg("e1"), pybind11::arg("e2"), pybind11::arg("p1"), pybind11::arg("p2"),
      pybind11::arg("van_der_waals_bond") = false,
      R"delim(
      Checks whether a bond exists between two atoms based on their distance.

      :param e1: Element type of first atom
      :param e2: Element type of second atom
      :param p1: Position of first atom
      :param p2: Position of second atom
    )delim");
  bond_detector.def_static("bond_exists",
                           pybind11::overload_cast<const ElementType&, const ElementType&, const Position&, const Position&,
                                                   const PeriodicBoundaries&, bool>(&BondDetector::bondExists),
                           pybind11::arg("e1"), pybind11::arg("e2"), pybind11::arg("p1"), pybind11::arg("p2"),
                           pybind11::arg("pbc"), pybind11::arg("van_der_waals_bond") = false,
                           R"delim(
      Checks whether a bond exists between two atoms based on their distance.

      :param e1: Element type of first atom
      :param e2: Element type of second atom
      :param p1: Position of first atom
      :param p2: Position of second atom
      :param pbc: Periodic Boundaries
    )delim");

  pybind11::class_<SolidStateBondDetector> solid_state_bond_detector(m, "SolidStateBondDetector",
                                                                     R"delim(
      A class to detect bonds based on interatomic distances.

      This detector can handle both pure solid state structures and heterogeneous systems.
      For the given solid state atoms, nearest neighbor bond orders are used with a margin of 0.1 Angstrom. If it is a
      heterogeneous system the conventional BondDetector (see 'BondDetector') is used to determine any bond involving at
      least one non-solid state atom.
      A binary decision on whether a bond exists (resulting in a bond order of 1.0) or not (yielding a bond order
      of 0.0) is made.
      Bonds across periodic boundaries can optionally be given the order of -1.0.
      Bonds between solid state atoms can optionally be determined with van der Waals bonds based on an optional flag.
      The van der Waals radii are specified in scine_utilities.ElementInfo
    )delim");
  solid_state_bond_detector.def_static("detect_bonds",
                                       pybind11::overload_cast<const AtomCollection&, const std::unordered_set<unsigned>&, bool>(
                                           &SolidStateBondDetector::detectBonds),
                                       pybind11::arg("atom_collection"), pybind11::arg("solid_state_indices"),
                                       pybind11::arg("solid_state_van_der_waals_bond") = false, R"delim(
      Generates a BondOrderCollection from an AtomCollection based on interatomic distances.)delim");
  solid_state_bond_detector.def_static(
      "detect_bonds", pybind11::overload_cast<const PeriodicSystem&, bool, bool>(&SolidStateBondDetector::detectBonds),
      pybind11::arg("periodic_system"), pybind11::arg("bonds_across_boundaries_negative") = false,
      pybind11::arg("solid_state_van_der_waals_bond") = false, R"delim(
      Generates a BondOrderCollection from a PeriodicSystem based on interatomic distances and periodic boundary conditions.)delim");
  solid_state_bond_detector.def_static(
      "detect_bonds",
      pybind11::overload_cast<const AtomCollection&, const PeriodicBoundaries&, const std::unordered_set<unsigned>&, bool, bool>(
          &SolidStateBondDetector::detectBonds),
      pybind11::arg("atom_collection"), pybind11::arg("pbc"), pybind11::arg("solid_state_indices"),
      pybind11::arg("bonds_across_boundaries_negative") = false, pybind11::arg("solid_state_van_der_waals_bond") = false, R"delim(
      Generates a BondOrderCollection from an AtomCollection based on interatomic distances and periodic boundary conditions.)delim");
  solid_state_bond_detector.def_static(
      "detect_bonds",
      pybind11::overload_cast<const ElementTypeCollection&, const PositionCollection&, const std::unordered_set<unsigned>&, bool>(
          &SolidStateBondDetector::detectBonds),
      pybind11::arg("elements"), pybind11::arg("positions"), pybind11::arg("solid_state_indices"),
      pybind11::arg("solid_state_van_der_waals_bond") = false, R"delim(
      Generates a BondOrderCollection from an ElementTypeCollection and a PositionCollection based on interatomic distances.)delim");
  solid_state_bond_detector.def_static(
      "detect_bonds",
      pybind11::overload_cast<const ElementTypeCollection&, const PositionCollection&, const PeriodicBoundaries&,
                              const std::unordered_set<unsigned>&, bool, bool>(&SolidStateBondDetector::detectBonds),
      pybind11::arg("elements"), pybind11::arg("positions"), pybind11::arg("pbc"), pybind11::arg("solid_state_indices"),
      pybind11::arg("bonds_across_boundaries_negative") = false, pybind11::arg("solid_state_van_der_waals_bond") = false, R"delim(
      Generates a BondOrderCollection from an ElementTypeCollection and a PositionCollection based on interatomic distances and periodic boundary conditions.)delim");
}
