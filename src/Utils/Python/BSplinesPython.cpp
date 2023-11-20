/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Math/BSplines/ReactionProfileInterpolation.h>
#include <Utils/Pybind.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine;

void init_bspline_functionalities(pybind11::module& m) {
  auto bsplines = m.def_submodule("bsplines");
  bsplines.doc() = R"delim(
    A collection of functions to generate and work with B-Splines of molecular
    trajectories.
  )delim";

  pybind11::class_<Utils::BSplines::TrajectorySpline> bspline{bsplines, "TrajectorySpline", R"delim(
      A class representing a B-Spline fit.
    )delim"};

  bspline.def(pybind11::init<Utils::ElementTypeCollection, Eigen::VectorXd, Eigen::MatrixXd>());

  bspline.def(pybind11::init<Utils::ElementTypeCollection, Eigen::VectorXd, Eigen::MatrixXd, double>());

  bspline.def("evaluate", &Utils::BSplines::TrajectorySpline::evaluate, pybind11::arg("position"), pybind11::arg("degree") = 0,
              R"delim(
      Fit a spline to all currently stored data.

      :param position: The position at which to evaluate the spline. has to be
                       in the interval [0,1]
      :param degree: The degree of the polynomial fit.
      :return: A Tuple of the fitted data: (energy: double , structure: AtomCollection).
    )delim");

  bspline.def_readonly("elements", &Utils::BSplines::TrajectorySpline::elements,
                       R"delim(
      The elements of the atoms that are represented.
    )delim");

  bspline.def_readonly("data", &Utils::BSplines::TrajectorySpline::data,
                       R"delim(
      The spline data to be fitted to.
    )delim");

  bspline.def_readonly("ts_position", &Utils::BSplines::TrajectorySpline::tsPosition,
                       R"delim(
      The TS position in the spline [0.0, 1.0], -1.0 if none is present.
    )delim");

  bspline.def_readonly("knots", &Utils::BSplines::TrajectorySpline::knots,
                       R"delim(
      The knot positions of the spline.
    )delim");

  pybind11::class_<Utils::BSplines::ReactionProfileInterpolation> rpi{bsplines, "ReactionProfileInterpolation",
                                                                      R"delim(
      Factory for B-Splines of reaction paths, including interpolation of an
      energy associated with the molecular structures.
    )delim"};

  rpi.def(pybind11::init<>());

  rpi.def("append_structure", &Utils::BSplines::ReactionProfileInterpolation::appendStructure, pybind11::arg("atoms"),
          pybind11::arg("energy"), pybind11::arg("is_the_transition_state") = false,
          R"delim(
      Adds a datapoint to the internal storage.
      Given atoms have to match those previously given, if there were any.

      :param atoms: An atom collection
      :param energy: The energy of the given atomic configuration.
    )delim");

  rpi.def("clear", &Utils::BSplines::ReactionProfileInterpolation::clear,
          R"delim(
      Clear all previously stored data.
    )delim");

  rpi.def("spline", &Utils::BSplines::ReactionProfileInterpolation::spline,
          pybind11::arg("n_interpolation_points") = 11, pybind11::arg("degree") = 3,
          R"delim(
      Fit a spline to all currently stored data.

      :param n_interpolation_points: Number of points to use for interpolation.
      :param degree: The maximum degree of the fit polynominals.
      :return: A BSpline object fitted to all stored data.
    )delim");

  rpi.def("current_ts_position", &Utils::BSplines::ReactionProfileInterpolation::getCurrentTSPosition,
          R"delim(
      :return: The current position of the transiton state in the intverval [0.0, 1.0]
    )delim");

  rpi.def_readwrite("use_quaternion_fit", &Utils::BSplines::ReactionProfileInterpolation::useQuaternionFit,
                    R"delim(
      If true, will determine distance along path by fitting each structure
      to the previous one. This fit is done in mass-weighted coordinates using
      quaternions.
    )delim");
}
