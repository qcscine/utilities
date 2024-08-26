/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/IO/TurbomoleMinimalBasisfile.h"
#include "Utils/Pybind.h"
#include "boost/exception/diagnostic_information.hpp"
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

void init_io(pybind11::module& m) {
  pybind11::module io = m.def_submodule("io");
  io.doc() = R"(The ``io`` submodule defines file input and output methods.)";

  io.def("read", pybind11::overload_cast<const std::string&>(&ChemicalFileHandler::read), "Reads an AtomCollection from a file");

  io.def("write",
         pybind11::overload_cast<const std::string&, const AtomCollection&, const std::string&>(&ChemicalFileHandler::write),
         pybind11::arg("filename"), pybind11::arg("atoms"), pybind11::arg("comment") = "",
         "Write an AtomCollection to a file");

  io.def("write_topology",
         pybind11::overload_cast<const std::string&, const AtomCollection&, const BondOrderCollection&, const std::string&>(
             &ChemicalFileHandler::write),
         pybind11::arg("filename"), pybind11::arg("atoms"), pybind11::arg("bondorders"), pybind11::arg("comment") = "",
         "Write an AtomCollection and a BondOrderCollection to a file");

  pybind11::enum_<MolecularTrajectoryIO::format>(io, "TrajectoryFormat")
      .value("Xyz", MolecularTrajectoryIO::format::xyz)
      .value("Binary", MolecularTrajectoryIO::format::binary)
      .value("Pdb", MolecularTrajectoryIO::format::pdb);

  io.def("read_trajectory",
         pybind11::overload_cast<const enum MolecularTrajectoryIO::format, const std::string&>(&MolecularTrajectoryIO::read),
         "Reads trajectory from a file");

  io.def("write_trajectory",
         pybind11::overload_cast<const enum MolecularTrajectoryIO::format, const std::string&, const MolecularTrajectory&>(
             &MolecularTrajectoryIO::write),
         pybind11::arg("f"), pybind11::arg("fileName"), pybind11::arg("m"), "Writes trajectory to a file");

  io.def("write_trajectory",
         pybind11::overload_cast<const enum MolecularTrajectoryIO::format, const std::string&,
                                 const MolecularTrajectory&, const BondOrderCollection&>(&MolecularTrajectoryIO::write),
         pybind11::arg("f"), pybind11::arg("fileName"), pybind11::arg("m"), pybind11::arg("bondorders"),
         "Writes trajectory to a file");

  io.def("read_trajectory",
         pybind11::overload_cast<const enum MolecularTrajectoryIO::format, const std::string&>(&MolecularTrajectoryIO::read),
         "Reads trajectory to a file");

  io.def("read_basis", &readTurbomoleBasisfile, "Reads basis file in Turbomole Format");
}
