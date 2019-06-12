/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/Atom.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

using namespace Scine::Utils;

void init_atom(pybind11::module& m) {
  pybind11::class_<Atom> atom(m, "Atom");
  atom.def(pybind11::init<ElementType, Position>(), pybind11::arg("e") = ElementType::none,
           pybind11::arg("p") = Position(0, 0, 0));
  atom.def_property("element", &Atom::getElementType, &Atom::setElementType, "The element type of the atom");
  atom.def_property("position", &Atom::getPosition, &Atom::setPosition, "The position of the atom");
}
