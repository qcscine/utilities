/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
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
  molecular_trajectory.def(
      "__getitem__",
      [&](MolecularTrajectory& traj, int i) -> PositionCollection {
        if (i >= traj.size()) {
          throw std::out_of_range("Given index is out of range");
        }
        if (i < 0) {
          if (-1 * i > traj.size()) {
            throw std::out_of_range("Given index is out of range");
          }
          return traj[traj.size() + i];
        }
        return traj[i];
      },
      "Access a PositionCollection frame of the trajectory.");
  molecular_trajectory.def(
      "__getitem__",
      [&](MolecularTrajectory& traj, pybind11::slice slice) -> MolecularTrajectory {
        auto energies = traj.getEnergies();
        bool addEnergies = energies.size() > 0;
        std::size_t start = 0, stop = 0, step = 0, slicelength = 0;
        if (!slice.compute(traj.size(), &start, &stop, &step, &slicelength)) {
          throw pybind11::error_already_set();
        }
        MolecularTrajectory sliced = MolecularTrajectory();
        sliced.setElementTypes(traj.getElementTypes());
        for (std::size_t i = 0; i < slicelength; ++i, start += step) {
          if (addEnergies) {
            sliced.push_back(traj.at(start), energies.at(start));
          }
          else {
            sliced.push_back(traj.at(start));
          }
        }
        return sliced;
      },
      "Access a sub-trajectory of the trajectory based on slicing.");
  molecular_trajectory.def(
      "__delitem__",
      [&](MolecularTrajectory& traj, int i) {
        const int N = traj.size();
        if (i >= N || -1 * i > N) {
          throw std::out_of_range("Cannot delete out of range index");
        }
        if (i < 0) {
          i = N + i;
        }
        auto energies = traj.getEnergies();
        bool addEnergies = energies.size() > 0;
        MolecularTrajectory result = MolecularTrajectory();
        result.setElementTypes(traj.getElementTypes());
        for (int j = 0; j < N; ++j) {
          if (j == i) {
            continue;
          }
          if (addEnergies) {
            result.push_back(traj.at(j), energies.at(j));
          }
          else {
            result.push_back(traj.at(j));
          }
        }
        traj = result;
      },
      "Allow the python delete function based on index.");

  molecular_trajectory.def(
      "__iter__", [](const MolecularTrajectory& m) { return pybind11::make_iterator(m.begin(), m.end()); },
      pybind11::keep_alive<0, 1>() // Keep object alive while iterator exists
  );

  molecular_trajectory.def("__len__", &MolecularTrajectory::size);

  // Operators
  molecular_trajectory.def(pybind11::self *= double());
  molecular_trajectory.def(pybind11::self /= double());
  molecular_trajectory.def(pybind11::self * double());
  molecular_trajectory.def(pybind11::self / double());
}
