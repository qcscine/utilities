/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/IO/TurbomoleMinimalBasisfile.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/MolecularTrajectory.h>
#include <pybind11/eigen.h>
#include <pybind11/eval.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

namespace {

void raiseDeprecationWarning() {
  pybind11::exec(
      R"delim(
      import warnings
      warnings.warn("'IO' submodule is deprecated in favor of PEP8 compliant 'io'", DeprecationWarning)
    )delim");
}

} // namespace

void init_io(pybind11::module& m) {
  pybind11::module io = m.def_submodule("io");
  io.doc() = R"(The ``io`` submodule defines file input and output methods.)";

  io.def("read", pybind11::overload_cast<const std::string&>(&ChemicalFileHandler::read), "Reads an AtomCollection from a file");

  io.def("write", pybind11::overload_cast<const std::string&, const AtomCollection&>(&ChemicalFileHandler::write),
         "Write an AtomCollection to a file");

  io.def("write_topology",
         pybind11::overload_cast<const std::string&, const AtomCollection&, const BondOrderCollection&>(&ChemicalFileHandler::write),
         "Write an AtomCollection and a BondOrderCollection to a file");

  pybind11::enum_<MolecularTrajectoryIO::format>(io, "TrajectoryFormat")
      .value("Xyz", MolecularTrajectoryIO::format::xyz)
      .value("Binary", MolecularTrajectoryIO::format::binary);

  io.def("write_trajectory",
         pybind11::overload_cast<const enum MolecularTrajectoryIO::format, const std::string&, const MolecularTrajectory&>(
             &MolecularTrajectoryIO::write),
         "Writes trajectory to a file");

  io.def("read_basis", &readTurbomoleBasisfile, "Reads basis file in Turbomole Format");

  /* Add deprecated IO module */
  pybind11::class_<ChemicalFileHandler> IO(m, "IO");
  IO.doc() = R"(Deprecated. Prefer lowercase io)";
  IO.def(
      "read",
      [](const std::string& filename) {
        raiseDeprecationWarning();
        return ChemicalFileHandler::read(filename);
      },
      "Reads an AtomCollection from a file");

  IO.def(
      "write",
      [](const std::string& a, const AtomCollection& b) {
        raiseDeprecationWarning();
        return ChemicalFileHandler::write(a, b);
      },
      "Write an AtomCollection to a file");

  io.def(
      "write_topology",
      [](const std::string& a, const AtomCollection& b, const BondOrderCollection& c) {
        raiseDeprecationWarning();
        return ChemicalFileHandler::write(a, b, c);
      },
      "Write an AtomCollection and a BondOrderCollection to a file");

  IO.def(
      "write_trajectory",
      [](const MolecularTrajectoryIO::format a, const std::string& b, const MolecularTrajectory& c) {
        raiseDeprecationWarning();
        return MolecularTrajectoryIO::write(a, b, c);
      },
      "Writes trajectory to a file");

  IO.def(
      "read_basis",
      [](const std::string& x) {
        raiseDeprecationWarning();
        return readTurbomoleBasisfile(x);
      },
      "Reads basis file in Turbomole Format");

  IO.attr("format") = io.attr("TrajectoryFormat");
  IO.attr("xyz") = io.attr("TrajectoryFormat").attr("Xyz");
  IO.attr("binary") = io.attr("TrajectoryFormat").attr("Binary");
}
