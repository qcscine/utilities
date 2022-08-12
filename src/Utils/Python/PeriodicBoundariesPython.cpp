/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/DataStructures/PeriodicBoundaries.h>
#include <Utils/Geometry/GeometryUtilities.h>
#include <Utils/Pybind.h>
#include <pybind11/operators.h>

using namespace Scine::Utils;

void init_periodic_boundaries(pybind11::module& m) {
  pybind11::class_<PeriodicBoundaries> periodic_boundaries(m, "PeriodicBoundaries", R"delim(
   Class to represent periodic boundary conditions and apply periodic operations on PositionCollections.
   The class can be initialized either directly via a matrix defining the imaged cell or via the lengths and angles of the unit cell.
   The cell matrix is defined with following properties:
   vector a is defined to be along the x-axis
   vector b is defined to be in the x-y-plane
   The transformation matrix is defined as:

   .. math::
      \begin{split}
     matrix &= \begin{bmatrix}
              a0 & a1 & a2 \\
              b0 & b1 & b2 \\
              c0 & c1 & c2
              \end{bmatrix} \\

    \alpha &= \sphericalangle (\vec{b}, \vec{c})\\

    \beta &= \sphericalangle (\vec{a}, \vec{c})\\

    \gamma &= \sphericalangle (\vec{a}, \vec{b})\\
      \end{split})delim");
  periodic_boundaries.def(pybind11::init<double, const std::string&>(), pybind11::arg("cubeLength") = 1.0,
                          pybind11::arg("periodicity") = "xyz",
                          "Initialize cubic periodic boundaries with the given side length.");
  periodic_boundaries.def(pybind11::init<const Eigen::Matrix3d, const std::string&>(), pybind11::arg("matrix"),
                          pybind11::arg("periodicity") = "xyz",
                          "Initialize periodic boundaries with a particular cell matrix.");
  periodic_boundaries.def(pybind11::init<const Eigen::Vector3d&, const Eigen::Vector3d&, bool, bool, const std::string&>(),
                          pybind11::arg("lengths"), pybind11::arg("angles"), pybind11::arg("isBohr") = true,
                          pybind11::arg("isDegrees") = true, pybind11::arg("periodicity") = "xyz",
                          "Initialize from lengths of cell vectors and angles between them");
  periodic_boundaries.def(
      pybind11::init<std::string, std::string, bool, bool>(), pybind11::arg("periodicBoundariesString"),
      pybind11::arg("delimiter") = ",", pybind11::arg("isBohr") = true, pybind11::arg("isDegrees") = true,
      "Initialize from lengths of cell vectors and angles between them that are written in a string "
      "and separated by some delimiter");

  periodic_boundaries.def_property_readonly("a", &PeriodicBoundaries::getA, "Unit cell vector a");
  periodic_boundaries.def_property_readonly("b", &PeriodicBoundaries::getB, "Unit cell vector b");
  periodic_boundaries.def_property_readonly("c", &PeriodicBoundaries::getC, "Unit cell vector c");

  periodic_boundaries.def_property_readonly("alpha", &PeriodicBoundaries::getAlpha, "Unit cell angle alpha between b and c");
  periodic_boundaries.def_property_readonly("beta", &PeriodicBoundaries::getBeta, "Unit cell angle beta between a and c");
  periodic_boundaries.def_property_readonly("gamma", &PeriodicBoundaries::getGamma, "Unit cell angle gamma between a and b");

  periodic_boundaries.def_property_readonly("lengths", &PeriodicBoundaries::getLengths, "Lengths of the three unit vectors");
  periodic_boundaries.def_property_readonly("angles", &PeriodicBoundaries::getAngles,
                                            "Angles between the three unit vectors in degrees");

  periodic_boundaries.def_property("matrix", &PeriodicBoundaries::getCellMatrix, &PeriodicBoundaries::setCellMatrix,
                                   "The underlying matrix governing the periodic boundaries.");

  periodic_boundaries.def_property("periodicity", &PeriodicBoundaries::getPeriodicity,
                                   pybind11::overload_cast<const std::vector<bool>>(&PeriodicBoundaries::setPeriodicity),
                                   "The periodicity of the cell.");

  periodic_boundaries.def(
      "is_ortho_rhombic", &PeriodicBoundaries::isOrthoRhombic, pybind11::arg("eps") = 1e-2,
      "Returns whether the cell is orthorhombic. The optional parameter gives the tolerance around 90 degrees.");

  periodic_boundaries.def(
      "transform", pybind11::overload_cast<const PositionCollection&, bool>(&PeriodicBoundaries::transform, pybind11::const_),
      pybind11::arg("positions"), pybind11::arg("relativeToCartesian") = true,
      "Get a transformed PositionCollection from Relative to Cartesian Coordinates if boolean set to True and vice "
      "versa if set to False.");

  periodic_boundaries.def("transform",
                          pybind11::overload_cast<const Position&, bool>(&PeriodicBoundaries::transform, pybind11::const_),
                          pybind11::arg("position"), pybind11::arg("relativeToCartesian") = true,
                          "Get a transformed Position from relative to Cartesian Coordinates if boolean set to True "
                          "and vice versa if set to False.");

  periodic_boundaries.def(
      "transform_in_place",
      pybind11::overload_cast<PositionCollection&, bool>(&PeriodicBoundaries::transformInPlace, pybind11::const_),
      pybind11::arg("positions"), pybind11::arg("relativeToCartesian") = true,
      "Transform given PositionCollection from relative to Cartesian Coordinates if boolean set to True and vice versa "
      "for set to False.");

  periodic_boundaries.def(
      "transform_in_place",
      pybind11::overload_cast<Eigen::Ref<Position>, bool>(&PeriodicBoundaries::transformInPlace, pybind11::const_),
      pybind11::arg("position"), pybind11::arg("relativeToCartesian") = true,
      "Transform given Position from relative to Cartesian Coordinates if boolean set to True and "
      "vice versa for set to False.");

  periodic_boundaries.def("translate_positions_into_cell",
                          pybind11::overload_cast<const PositionCollection&, const Eigen::RowVector3d&>(
                              &PeriodicBoundaries::translatePositionsIntoCell, pybind11::const_),
                          pybind11::arg("positions"), Arg("relShift") = Position(0, 0, 0),
                          "Get a PositionCollection translated into the unit cell. Optionally you can give an "
                          "additional shift vector in Relative Coordinates.");

  periodic_boundaries.def("translate_positions_into_cell",
                          pybind11::overload_cast<const Position&, const Eigen::RowVector3d&>(
                              &PeriodicBoundaries::translatePositionsIntoCell, pybind11::const_),
                          pybind11::arg("positions"), Arg("relShift") = Position(0, 0, 0),
                          "Get a Position translated into the unit cell. Optionally you can give an additional shift "
                          "vector in Relative Coordinates.");

  periodic_boundaries.def("translate_positions_into_cell_in_place",
                          pybind11::overload_cast<Eigen::Ref<Position>, const Eigen::RowVector3d&>(
                              &PeriodicBoundaries::translatePositionsIntoCellInPlace, pybind11::const_),
                          pybind11::arg("positions"), Arg("relShift") = Position(0, 0, 0),
                          "Translate given Position into the unit cell. Optionally you can give an additional shift "
                          "vector in Relative Coordinates.");

  periodic_boundaries.def("translate_positions_into_cell_in_place",
                          pybind11::overload_cast<PositionCollection&, const Eigen::RowVector3d&>(
                              &PeriodicBoundaries::translatePositionsIntoCellInPlace, pybind11::const_),
                          pybind11::arg("positions"), Arg("relShift") = Position(0, 0, 0),
                          "Translate given PositionCollection into the unit cell. Optionally you can give an "
                          "additional shift vector in Relative Coordinates.");

  periodic_boundaries.def("__str__", &PeriodicBoundaries::getPeriodicBoundariesString, pybind11::arg("delimiter") = ",",
                          "String of all cell lengths and angles.");

  // operators
  periodic_boundaries.def(pybind11::self == pybind11::self);
  periodic_boundaries.def(pybind11::self != pybind11::self);
  periodic_boundaries.def(
      "__mul__", [&](const PeriodicBoundaries& pbc, double x) { return pbc * x; }, pybind11::is_operator());
  periodic_boundaries.def(
      "__imul__", [&](PeriodicBoundaries& pbc, double x) { return pbc *= x; }, pybind11::is_operator());
  periodic_boundaries.def(
      "__mul__", [&](const PeriodicBoundaries& pbc, const Eigen::Vector3d& x) { return pbc * x; }, pybind11::is_operator());
  periodic_boundaries.def(
      "__imul__", [&](PeriodicBoundaries& pbc, const Eigen::Vector3d& x) { return pbc *= x; }, pybind11::is_operator());
}
