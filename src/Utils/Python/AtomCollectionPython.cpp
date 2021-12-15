/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/AtomCollection.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_atom_collection(pybind11::module& m) {
  pybind11::class_<AtomCollection> atom_collection(m, "AtomCollection");
  atom_collection.def(pybind11::init<int>(), pybind11::arg("N") = 0, "Initialize a particular number of empty atoms");
  atom_collection.def(pybind11::init<ElementTypeCollection, PositionCollection>(),
                      "Initialize from element types and positions");

  // OO pairs
  atom_collection.def_property("elements", &AtomCollection::getElements, &AtomCollection::setElements,
                               "All element types of the collection");

  atom_collection.def("set_element", &AtomCollection::setElement, "Set an element type");
  atom_collection.def("get_element", &AtomCollection::getElement, "Get an element type");

  atom_collection.def_property("positions", &AtomCollection::getPositions, &AtomCollection::setPositions,
                               "All positions of the collection");

  atom_collection.def("set_position", &AtomCollection::setPosition, "Set a position");
  atom_collection.def("get_position", &AtomCollection::getPosition, "Get a position");

  // std::vector-like functions
  atom_collection.def("size", &AtomCollection::size, "Get how many atoms are in the collection");
  atom_collection.def("clear", &AtomCollection::clear, "Remove all atoms from the collection");
  atom_collection.def("resize", &AtomCollection::resize, "Resize the collection. Does not preserve contained data.");
  atom_collection.def("at", &AtomCollection::at, "Access a particular atom");
  atom_collection.def("push_back", &AtomCollection::push_back, "Add a new atom");

  // Comparison operators
  atom_collection.def(pybind11::self == pybind11::self);
  atom_collection.def(pybind11::self != pybind11::self);

  // Addition operators
  atom_collection.def(pybind11::self + pybind11::self);
  atom_collection.def(pybind11::self += pybind11::self);

  // Sequence magic methods
  atom_collection.def(
      "__iter__", [](const AtomCollection& s) { return pybind11::make_iterator(s.begin(), s.end()); },
      pybind11::keep_alive<0, 1>() // Keep object alive while iterator exists
  );
  atom_collection.def("__len__", &AtomCollection::size);
  atom_collection.def("__getitem__", &AtomCollection::operator[]);
  atom_collection.def("__getitem__", [&](AtomCollection& coll, pybind11::slice slice) -> AtomCollection {
    std::size_t start = 0, stop = 0, step = 0, slicelength = 0;
    if (!slice.compute(coll.size(), &start, &stop, &step, &slicelength)) {
      throw pybind11::error_already_set();
    }
    AtomCollection sliced(slicelength);
    for (std::size_t i = 0; i < slicelength; ++i, start += step) {
      sliced.setElement(i, coll.getElement(start));
      sliced.setPosition(i, coll.getPosition(start));
    }
    return sliced;
  });
  atom_collection.def("__delitem__", [&](AtomCollection& coll, int i) {
    const int N = coll.size();
    if (i >= N) {
      throw std::out_of_range("Cannot delete out of range index");
    }
    AtomCollection result(N - 1);
    for (int j = 0; j < i; ++j) {
      result.setElement(j, coll.getElement(j));
      result.setPosition(j, coll.getPosition(j));
    }
    for (int j = i; j < N - 1; ++j) {
      result.setElement(j, coll.getElement(j + 1));
      result.setPosition(j, coll.getPosition(j + 1));
    }
    coll = result;
  });

  // Copy support
  atom_collection.def("__copy__", [](const AtomCollection& c) -> AtomCollection { return AtomCollection(c); });
  atom_collection.def("__deepcopy__", [](const AtomCollection& c, pybind11::dict /* memo */) -> AtomCollection {
    return AtomCollection(c);
  });
}
