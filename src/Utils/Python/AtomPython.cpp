/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/Atom.h>
#include <Utils/Pybind.h>
#include <pybind11/stl_bind.h>

using namespace Scine::Utils;

void init_atom(pybind11::module& m) {
  pybind11::class_<Atom> atom(m, "Atom");
  atom.def(pybind11::init<ElementType, Position>(), Arg("e") = ElementType::none, Arg("p") = Position(0, 0, 0));
  atom.def_property("element", &Atom::getElementType, &Atom::setElementType, "The element type of the atom");
  atom.def_property("position", &Atom::getPosition, &Atom::setPosition, "The position of the atom");
}
