/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/MolecularTrajectory.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_molecular_trajectory(pybind11::module& m) {
  pybind11::class_<MolecularTrajectory> molecular_trajectory(m, "MolecularTrajectory");

  molecular_trajectory.def(pybind11::init<>());

  molecular_trajectory.def_property("elements", &MolecularTrajectory::getElementTypes,
                                    &MolecularTrajectory::setElementTypes, "Element types of the atoms");

  molecular_trajectory.def("set_element_type", &MolecularTrajectory::setElementType, "Set a single element type");

  molecular_trajectory.def("clear", &MolecularTrajectory::clear, "Removes all steps in the trajectory, but not element types");

  molecular_trajectory.def("clear_energies", &MolecularTrajectory::clearEnergies, "Removes all energies of the trajectory");

  molecular_trajectory.def("resize", &MolecularTrajectory::resize, "Sets the number of PositionCollections in the trajectory");

  molecular_trajectory.def("push_back", pybind11::overload_cast<PositionCollection>(&MolecularTrajectory::push_back),
                           pybind11::arg("positions"), "Add a new set of positions to the trajectory");

  molecular_trajectory.def("push_back", pybind11::overload_cast<PositionCollection, double>(&MolecularTrajectory::push_back),
                           pybind11::arg("positions"), pybind11::arg("energy"),
                           "Add a new set of positions to the trajectory with its corresponding energy");

  molecular_trajectory.def("empty", &MolecularTrajectory::empty, "Returns whether no structures are present");

  molecular_trajectory.def("size", &MolecularTrajectory::size, "Returns the number of structures in the trajectory");

  molecular_trajectory.def("molecular_size", &MolecularTrajectory::molecularSize, "Returns the number of atoms in the structure");

  molecular_trajectory.def("get_energies", &MolecularTrajectory::getEnergies, "Returns the energies of the trajectory");

  // Sequence magic methods
  molecular_trajectory.def("__getitem__", pybind11::overload_cast<int>(&MolecularTrajectory::operator[], pybind11::const_),
                           "Access a PositionCollection frame of the trajectory");

  molecular_trajectory.def("__iter__",
                           [](const MolecularTrajectory& m) { return pybind11::make_iterator(m.begin(), m.end()); },
                           pybind11::keep_alive<0, 1>() // Keep object alive while iterator exists
  );

  molecular_trajectory.def("__len__", &MolecularTrajectory::size);

  // Operators
  molecular_trajectory.def(pybind11::self *= double());
  molecular_trajectory.def(pybind11::self /= double());
  molecular_trajectory.def(pybind11::self * double());
  molecular_trajectory.def(pybind11::self / double());
}
