/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_io(pybind11::module& m) {
  pybind11::class_<ChemicalFileHandler> file_handler(m, "IO");

  file_handler.def_static("read", pybind11::overload_cast<const std::string&>(&ChemicalFileHandler::read),
                          "Reads an AtomCollection from a file");

  file_handler.def_static("write",
                          pybind11::overload_cast<const std::string&, const AtomCollection&>(&ChemicalFileHandler::write),
                          "Write an AtomCollection to a file");

  file_handler.def_static("write_topology",
                          pybind11::overload_cast<const std::string&, const AtomCollection&, const BondOrderCollection&>(
                              &ChemicalFileHandler::write),
                          "Write an AtomCollection and a BondOrderCollection to a file");
}
