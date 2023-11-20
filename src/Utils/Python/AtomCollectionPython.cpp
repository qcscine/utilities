/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>
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

  atom_collection.def_property("residues", &AtomCollection::getResidues, &AtomCollection::setResidues,
                               "The residue information (residue label, chain label, residue index)");

  atom_collection.def("set_residue_info", &AtomCollection::setResidueInformation,
                      "Set the residue information (residue label, chain label, residue index)");
  atom_collection.def("get_residue_info", &AtomCollection::getResidueInformation,
                      "Get the residue information (residue label, chain label, residue index)");

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
  atom_collection.def("append", &AtomCollection::push_back, "Add a new atom");
  atom_collection.def("swap_indices", &AtomCollection::swapIndices, pybind11::arg("i"), pybind11::arg("j"),
                      "Swap two atoms by the specified indices");

  // Comparison operators
  atom_collection.def(pybind11::self == pybind11::self);
  atom_collection.def(pybind11::self != pybind11::self);
  atom_collection.def("is_approx", &AtomCollection::isApprox, pybind11::arg("other_system"),
                      pybind11::arg("epsilon") = 1e-6, "Allows to set the accuracy of the fuzzy comparison.");
  atom_collection.def(
      "__lt__", [&](const AtomCollection& lhs, const AtomCollection& rhs) -> bool { return lhs.size() < rhs.size(); },
      pybind11::arg("other_system"), "AtomCollections are compared based on their size");
  atom_collection.def(
      "__le__", [&](const AtomCollection& lhs, const AtomCollection& rhs) -> bool { return lhs.size() <= rhs.size(); },
      pybind11::arg("other_system"), "AtomCollections are compared based on their size");
  atom_collection.def(
      "__gt__", [&](const AtomCollection& lhs, const AtomCollection& rhs) -> bool { return lhs.size() > rhs.size(); },
      pybind11::arg("other_system"), "AtomCollections are compared based on their size");
  atom_collection.def(
      "__ge__", [&](const AtomCollection& lhs, const AtomCollection& rhs) -> bool { return lhs.size() >= rhs.size(); },
      pybind11::arg("other_system"), "AtomCollections are compared based on their size");

  // Addition operators
  atom_collection.def(pybind11::self + pybind11::self);
  atom_collection.def(pybind11::self += pybind11::self);

  // Sequence magic methods
  atom_collection.def(
      "__iter__", [](const AtomCollection& s) { return pybind11::make_iterator(s.begin(), s.end()); },
      pybind11::keep_alive<0, 1>() // Keep object alive while iterator exists
  );
  atom_collection.def("__len__", &AtomCollection::size);
  atom_collection.def(
      "__getitem__",
      [&](AtomCollection& coll, int i) -> Atom {
        if (i >= coll.size()) {
          throw std::out_of_range("Given index is out of range");
        }
        if (i < 0) {
          if (-1 * i > coll.size()) {
            throw std::out_of_range("Given index is out of range");
          }
          return coll[coll.size() + i];
        }
        return coll[i];
      },
      "Access an Atom of the AtomCollection.");
  atom_collection.def(
      "__getitem__",
      [&](AtomCollection& coll, pybind11::slice slice) -> AtomCollection {
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
      },
      "Access a sub-AtomCollection of the AtomCollection based on slicing.");
  atom_collection.def(
      "__delitem__",
      [&](AtomCollection& coll, int i) {
        const int N = coll.size();
        if (i >= N || -1 * i > N) {
          throw std::out_of_range("Cannot delete out of range index");
        }
        if (i < 0) {
          i = N + i;
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
      },
      "Allow the python delete function based on index.");

  // Copy support
  atom_collection.def("__copy__", [](const AtomCollection& c) -> AtomCollection { return AtomCollection(c); });
  atom_collection.def("__deepcopy__", [](const AtomCollection& c, pybind11::dict /* memo */) -> AtomCollection {
    return AtomCollection(c);
  });
  atom_collection.def(pybind11::pickle(
      [](const AtomCollection& ac) { // __getstate__
        /* Return a tuple that fully encodes the state of the object */
        std::vector<std::string> elements;
        for (const auto& e : ac.getElements()) {
          elements.push_back(ElementInfo::symbol(e));
        }
        return pybind11::make_tuple(elements, ac.getPositions());
      },
      [](pybind11::tuple t) { // __setstate__
        if (t.size() != 2)
          throw std::runtime_error("Invalid state for AtomCollection!");
        ElementTypeCollection elements;
        for (const auto& e : t[0].cast<std::vector<std::string>>()) {
          elements.push_back(ElementInfo::elementTypeForSymbol(e));
        }
        AtomCollection ac(elements, t[1].cast<PositionCollection>());
        return ac;
      }));
}
