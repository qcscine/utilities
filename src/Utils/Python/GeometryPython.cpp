/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;
using namespace Scine::Utils::Geometry;

void init_geometry(pybind11::module& m) {
  auto geometry_submodule = m.def_submodule("geometry");

  // This is going to fail, how to overload this?
  geometry_submodule.def(
      "translate_positions",
      pybind11::overload_cast<const PositionCollection&, const Eigen::Ref<Eigen::RowVector3d>&>(&translatePositions),
      "Translates a set of positions by a Displacement");

  geometry_submodule.def(
      "rotate_positions",
      pybind11::overload_cast<const PositionCollection&, const Eigen::Ref<Eigen::RowVector3d>&,
                              const Eigen::Ref<Eigen::RowVector3d>&, const Eigen::Ref<Eigen::RowVector3d>&>(&rotatePositions),
      "Rotate a set of positions according to two given vectors");

  geometry_submodule.def("rotate_positions",
                         pybind11::overload_cast<const PositionCollection&, const Eigen::Ref<Eigen::RowVector3d>&,
                                                 double, const Eigen::Ref<Eigen::RowVector3d>&>(&rotatePositions),
                         "Rotate a set of positions around a given axis");

  geometry_submodule.def("random_displacement",
                         pybind11::overload_cast<const PositionCollection&, const double, const double>(&randomDisplacement),
                         "Randomly displace a set of positions with a maximum displacement per coordinate. The "
                         "positions and the maximum displacement have to be given as arguments. The default seed is "
                         "42.0 and can be given as the third argument if needed.",
                         pybind11::arg("positions"), pybind11::arg("maxDisplacement"), pybind11::arg("seed") = 42.0);

  geometry_submodule.def("random_displacement_trajectory",
                         pybind11::overload_cast<const AtomCollection&, const int, const double>(&randomDisplacementTrajectory),
                         "Randomly displaces positions in AtomCollection and saves them as MolecularTrajectory. The "
                         "AtomCollection, the number of frames in the resulting trajectory and the maximum "
                         "displacement for each coordinate have to be given as arguments.",
                         pybind11::arg("atoms"), pybind11::arg("numFrames"), pybind11::arg("maxDisplacement"));

  geometry_submodule.def("get_index_of_closest_atom", &getIndexOfClosestAtom,
                         "Get the index of an atom that is closes to a spatial position");

  geometry_submodule.def("position_vector_to_matrix", &positionVectorToMatrix,
                         "Transform a 3N dimensional vector into a Nx3 matrix");

  geometry_submodule.def("position_matrix_to_vector", &positionMatrixToVector, "Transform a Nx3 matrix into a 3N vector");

  geometry_submodule.def("align_positions",
                         pybind11::overload_cast<const PositionCollection&, PositionCollection&>(&alignPositions),
                         "Rotate and translate positions to match a reference as closely as possible");
  geometry_submodule.def(
      "align_positions",
      pybind11::overload_cast<const PositionCollection&, PositionCollection&, const ElementTypeCollection&>(&alignPositions),
      "Rotate and translate positions to match a reference as closely as possible, uses masses as weight for the fit.");

  geometry_submodule.def("get_masses", &getMasses, "Get a vector of all element type masses");

  geometry_submodule.def("get_center_of_mass", pybind11::overload_cast<const AtomCollection&>(&getCenterOfMass),
                         "Get the center of mass position of an atom collection");

  geometry_submodule.def("get_center_of_mass",
                         pybind11::overload_cast<const PositionCollection&, const std::vector<double>&>(&getCenterOfMass),
                         "Get the center of mass position of a PositionCollection and a vector of masses");

  geometry_submodule.def("get_average_position", &getAveragePosition, "Calculate the average position of an atom collection");

  geometry_submodule.def("calculate_inertia_tensor", &calculateInertiaTensor, "Calculate the inertia tensor");

  pybind11::class_<PrincipalMomentsOfInertia> principal_moments_of_inertia(geometry_submodule,
                                                                           "PrincipalMomentsOfInertia");

  principal_moments_of_inertia.def_property_readonly("eigenvalues",
                                                     [](const PrincipalMomentsOfInertia& o) { return o.eigenvalues; });
  principal_moments_of_inertia.def_property_readonly("eigenvectors",
                                                     [](const PrincipalMomentsOfInertia& o) { return o.eigenvectors; });

  geometry_submodule.def("principal_inertial_moments", &calculatePrincipalMoments, "Calculate the principal moments of inertia");
}
