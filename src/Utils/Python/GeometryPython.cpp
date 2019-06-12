/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;
using namespace Scine::Utils::Geometry;

void init_geometry(pybind11::module& m) {
  auto geometry_submodule = m.def_submodule("Geometry");

  // This is going to fail, how to overload this?
  geometry_submodule.def("translate_positions",
                         pybind11::overload_cast<const PositionCollection&, const Displacement&>(&translatePositions),
                         "Translates a set of positions by a Displacement");

  geometry_submodule.def("get_index_of_closest_atom", &getIndexOfClosestAtom,
                         "Get the index of an atom that is closes to a spatial position");

  geometry_submodule.def("position_vector_to_matrix", &positionVectorToMatrix,
                         "Transform a 3N dimensional vector into a Nx3 matrix");

  geometry_submodule.def("position_matrix_to_vector", &positionMatrixToVector, "Transform a Nx3 matrix into a 3N vector");

  geometry_submodule.def("align_positions", &alignPositions,
                         "Rotate and translate positions to match a reference as closely as possible");

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
