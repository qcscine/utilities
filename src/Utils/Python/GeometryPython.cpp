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
using namespace Scine::Utils::Geometry::Distances;
using namespace Scine::Utils::Geometry::Manipulations;
using namespace Scine::Utils::Geometry::Properties;
using namespace Scine::Utils::Geometry::Transformations;

void init_geometry(pybind11::module& m) {
  auto geometry_submodule = m.def_submodule("geometry");

  geometry_submodule.def("translate_positions", &Manipulations::translatePositions, pybind11::arg("positions"),
                         pybind11::arg("translation"), "Translates a set of positions by a Displacement");

  geometry_submodule.def(
      "rotate_positions",
      pybind11::overload_cast<const PositionCollection&, const Eigen::RowVector3d&, const Eigen::RowVector3d&, const Eigen::RowVector3d&>(
          &Manipulations::rotatePositions),
      "Rotate a set of positions according to two given vectors");

  geometry_submodule.def(
      "rotate_positions",
      pybind11::overload_cast<const PositionCollection&, const Eigen::RowVector3d&, double, const Eigen::RowVector3d&>(
          &Manipulations::rotatePositions),
      "Rotate a set of positions around a given axis");

  geometry_submodule.def("random_displacement", &Manipulations::randomDisplacement,
                         "Randomly displace a set of positions with a maximum displacement per coordinate. The "
                         "positions and the maximum displacement have to be given as arguments.",
                         pybind11::arg("positions"), pybind11::arg("maxDisplacement"));

  geometry_submodule.def(
      "random_displacement_trajectory",
      pybind11::overload_cast<const AtomCollection&, unsigned int, double>(&Manipulations::randomDisplacementTrajectory),
      "Randomly displaces positions in AtomCollection and saves them as MolecularTrajectory. The "
      "AtomCollection, the number of frames in the resulting trajectory and the maximum "
      "displacement for each coordinate have to be given as arguments. Further, a specific seed can be specified. "
      "The default is 42.",
      pybind11::arg("atoms"), pybind11::arg("numFrames"), pybind11::arg("maxDisplacement"));

  geometry_submodule.def(
      "random_displacement_trajectory",
      pybind11::overload_cast<const AtomCollection&, unsigned int, double, unsigned int>(&Manipulations::randomDisplacementTrajectory),
      "Randomly displaces positions in AtomCollection and saves them as MolecularTrajectory. The "
      "AtomCollection, the number of frames in the resulting trajectory and the maximum "
      "displacement for each coordinate have to be given as arguments. Further, a specific seed can be specified to "
      "reproduce previous results.",
      pybind11::arg("atoms"), pybind11::arg("numFrames"), pybind11::arg("maxDisplacement"), pybind11::arg("seed"));

  geometry_submodule.def(
      "get_index_of_closest_atom",
      pybind11::overload_cast<const PositionCollection&, const Position&, double>(&Distances::getIndexOfClosestAtom),
      pybind11::arg("positions"), pybind11::arg("targetPosition"), pybind11::arg("squaredDistanceConsideredZero") = -1.0,
      "Get the index of an atom that is closest to a spatial position. ");

  geometry_submodule.def(
      "get_index_of_closest_atom",
      pybind11::overload_cast<const PositionCollection&, const Position&, const PeriodicBoundaries&, double>(
          &Distances::getIndexOfClosestAtom),
      pybind11::arg("positions"), pybind11::arg("targetPosition"), pybind11::arg("pbc"),
      pybind11::arg("squaredDistanceConsideredZero") = -1.0,
      "Get the index of an atom that is closest to a spatial position considering periodic boundary conditions. ");

  geometry_submodule.def("position_vector_to_matrix", &Transformations::positionVectorToMatrix,
                         "Transform a 3N dimensional vector into a Nx3 matrix");

  geometry_submodule.def("position_matrix_to_vector", &Transformations::positionMatrixToVector,
                         "Transform a Nx3 matrix into a 3N vector");

  geometry_submodule.def("align_positions",
                         pybind11::overload_cast<const PositionCollection&, PositionCollection&>(&Manipulations::alignPositions),
                         "Rotate and translate positions to match a reference as closely as possible");
  geometry_submodule.def(
      "align_positions",
      pybind11::overload_cast<const PositionCollection&, PositionCollection&, const ElementTypeCollection&>(
          &Manipulations::alignPositions),
      "Rotate and translate positions to match a reference as closely as possible, uses masses as weight for the fit.");

  geometry_submodule.def("get_masses", &Properties::getMasses, "Get a vector of all element type masses");

  geometry_submodule.def("get_center_of_mass", pybind11::overload_cast<const AtomCollection&>(&Properties::getCenterOfMass),
                         "Get the center of mass position of an atom collection");

  geometry_submodule.def(
      "get_center_of_mass",
      pybind11::overload_cast<const PositionCollection&, const std::vector<double>&>(&Properties::getCenterOfMass),
      "Get the center of mass position of a PositionCollection and a vector of masses");

  geometry_submodule.def("get_average_position", &Properties::getAveragePosition,
                         "Calculate the average position of an atom collection");

  geometry_submodule.def("calculate_inertia_tensor", &Properties::calculateInertiaTensor, "Calculate the inertia tensor");

  pybind11::class_<Properties::PrincipalMomentsOfInertia> principal_moments_of_inertia(geometry_submodule,
                                                                                       "PrincipalMomentsOfInertia");

  principal_moments_of_inertia.def_property_readonly(
      "eigenvalues", [](const Properties::PrincipalMomentsOfInertia& o) { return o.eigenvalues; });
  principal_moments_of_inertia.def_property_readonly(
      "eigenvectors", [](const Properties::PrincipalMomentsOfInertia& o) { return o.eigenvectors; });

  geometry_submodule.def("principal_inertial_moments", &Properties::calculatePrincipalMoments,
                         "Calculate the principal moments of inertia");

  geometry_submodule.def("distance",
                         pybind11::overload_cast<const PositionCollection&, const PositionCollection&>(&Distances::distance),
                         pybind11::arg("p1"), pybind11::arg("p2"),
                         R"delim(
      Get the distance between two PositionCollections.

      :param p1: First PositionCollection
      :param p2: Second PositionCollection
    )delim");

  geometry_submodule.def(
      "distance",
      pybind11::overload_cast<const PositionCollection&, const PositionCollection&, const PeriodicBoundaries&>(&Distances::distance),
      pybind11::arg("p1"), pybind11::arg("p2"), pybind11::arg("pbc"),
      R"delim(
      Get the distance between two PositionCollections with Periodic Boundaries considered.

      :param p1: First PositionCollection
      :param p2: Second PositionCollection
      :param pbc: Periodic Boundaries
    )delim");

  geometry_submodule.def("distance", pybind11::overload_cast<const Position&, const Position&>(&Distances::distance),
                         pybind11::arg("p1"), pybind11::arg("p2"),
                         R"delim(
      Get the distance between two Positions.

      :param p1: First Position
      :param p2: Second Position
    )delim");

  geometry_submodule.def(
      "distance", pybind11::overload_cast<const Position&, const Position&, const PeriodicBoundaries&>(&Distances::distance),
      pybind11::arg("p1"), pybind11::arg("p2"), pybind11::arg("pbc"),
      R"delim(
      Get the distance between two Positions with Periodic Boundaries considered.

      :param p1: First Position
      :param p2: Second Position
      :param pbc: Periodic Boundaries
    )delim");

  geometry_submodule.def("nearest_neighbors_bond_orders",
                         pybind11::overload_cast<const PositionCollection&, double>(&Distances::nearestNeighborsBondOrders),
                         pybind11::arg("positions"), pybind11::arg("margin") = 0.1,
                         R"delim(
      Get the BondOrderCollection of a Geometry where only nearest neighbors are connected by single bonds.

      The criterion for a bond is that at least ONE of the two positions has to be a nearest neighbor of the other
      one, and the nearest neighbor property does NOT have to be mutual. This means that the resulting bonding partner(s)
      of a position can be different to its nearest neighbors by the position being a nearest neighbor of its bonding
      partner(s).
      in short: if x is NN of y OR y is a NN of x: bond between x and y
      This also means that each position has at least one bond.

      :param positions: The PositionCollection
      :param margin: The margin to consider two distances to be equal
    )delim");

  geometry_submodule.def("nearest_neighbors_bond_orders",
                         pybind11::overload_cast<const PositionCollection&, const PeriodicBoundaries&, double>(
                             &Distances::nearestNeighborsBondOrders),
                         pybind11::arg("positions"), pybind11::arg("pbc"), pybind11::arg("margin") = 0.1,
                         R"delim(
      Get the BondOrderCollection of a Geometry where only nearest neighbors are connected by single bonds with Periodic Boundaries considered.

      The criterion for a bond is that at least ONE of the two positions has to be a nearest neighbor of the other
      one, and the nearest neighbor property does NOT have to be mutual. This means that the resulting bonding partner(s)
      of a position can be different to its nearest neighbors by the position being a nearest neighbor of its bonding
      partner(s).
      in short: if x is NN of y OR y is a NN of x: bond between x and y
      This also means that each position has at least one bond.

      :param positions: The PositionCollection
      :param pbc: Periodic Boundaries
      :param margin: The margin to consider two distances to be equal
    )delim");

  geometry_submodule.def("count_all_nearest_neighbors",
                         pybind11::overload_cast<const PositionCollection&, double>(&Distances::countAllNearestNeighbors),
                         pybind11::arg("positions"), pybind11::arg("margin") = 0.1,
                         R"delim(
      Count the number of the nearest neighbors for all Positions within a PositionCollection.

      :param positions: The PositionCollection
      :param margin: The margin to consider two distances to be equal
    )delim");

  geometry_submodule.def("count_all_nearest_neighbors",
                         pybind11::overload_cast<const PositionCollection&, const PeriodicBoundaries&, double>(
                             &Distances::countAllNearestNeighbors),
                         pybind11::arg("positions"), pybind11::arg("pbc"), pybind11::arg("margin") = 0.1,
                         R"delim(
      Count the number of the nearest neighbors for all Positions within a PositionCollection with Periodic Boundaries considered.

      :param positions: The PositionCollection
      :param pbc: Periodic Boundaries
      :param margin: The margin to consider two distances to be equal
    )delim");
}
