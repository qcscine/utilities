/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Bonds/BondDetector.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/Geometry/AtomCollection.h>
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

      The covalent radii were extracted from the Cambridge Structural Database (CSD) on 04/08/2020:
      https://www.ccdc.cam.ac.uk/support-and-resources/ccdcresources/Elemental_Radii.xlsx

      References:
      E. C. Meng, R. A. Lewis, Comput. Chem. 1991, 12, 891-898. [DOI: 10.1002/jcc.540120716]
      C. R. Groom, I. J. Bruno, M. P. Lightfoot and S. C. Ward, Acta Cryst. 2016, B72, 171-179. [DOI: 10.1107/S2052520616003954]
    )delim");
  bond_detector.def_static("detect_bonds", pybind11::overload_cast<const AtomCollection&>(&BondDetector::detectBonds),
                           pybind11::arg("atom_collection"), R"delim(
      Generates a BondOrderCollection from an AtomCollection based on interatomic distances.)delim");
  bond_detector.def_static(
      "detect_bonds",
      pybind11::overload_cast<const ElementTypeCollection&, const PositionCollection&>(&BondDetector::detectBonds),
      pybind11::arg("elements"), pybind11::arg("positions"), R"delim(
      Generates a BondOrderCollection from an ElementTypeCollection and a PositionCollection based on interatomic distances.)delim");
  bond_detector.def_static("bond_exists", &BondDetector::bondExists, pybind11::arg("e1"), pybind11::arg("e2"),
                           pybind11::arg("p1"), pybind11::arg("p2"),
                           R"delim(
      Checks whether a bond exists between two atoms based on their distance.

      :param e1: Element type of first atom
      :param e2: Element type of second atom
      :param p1: Position of first atom
      :param p2: Position of second atom
    )delim");
}
