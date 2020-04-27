/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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
  pybind11::class_<BondDetector> bond_detector(m, "BondDetector");
  bond_detector.def_static("detect_bonds", pybind11::overload_cast<const AtomCollection&>(&BondDetector::detectBonds), R"delim(
      Generates a BondOrderCollection from an AtomCollection based on covalent radii.
      A bond is detected if the distance between two atoms is smaller than the sum of their covalent radii plus 0.4 Angstrom.)delim");
  bond_detector.def_static(
      "detect_bonds",
      pybind11::overload_cast<const ElementTypeCollection&, const PositionCollection&>(&BondDetector::detectBonds), R"delim(
      Generates a BondOrderCollection from an ElementTypeCollection and a PositionCollection based on covalent radii.
      A bond is detected if the distance between two atoms is smaller than the sum of their covalent radii plus 0.4 Angstrom.)delim");
  bond_detector.def_static("bond_exists", &BondDetector::bondExists, R"delim(
      Checks whether a bond exists between two atoms based on covalent radii.
      A bond is detected if the distance between the two atoms is smaller than the sum of their covalent radii plus 0.4 Angstrom.)delim");
}