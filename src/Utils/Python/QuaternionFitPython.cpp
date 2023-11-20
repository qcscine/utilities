/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry.h>
#include <Utils/Math/QuaternionFit.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

using namespace Scine::Utils;

void init_quaternion_fit(pybind11::module& m) {
  pybind11::class_<QuaternionFit> quaternion_fit(m, "QuaternionFit");

  quaternion_fit.def(pybind11::init<Eigen::MatrixXd, Eigen::MatrixXd, bool>(), pybind11::arg("ref_matrix"),
                     pybind11::arg("fit_matrix"), pybind11::arg("improperRotationIsAllowed") = false,
                     "Initialize a QuaternionFit with a reference, "
                     "and a matrix to be fitted ");

  quaternion_fit.def(pybind11::init<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd, bool>(),
                     pybind11::arg("ref_matrix"), pybind11::arg("fit_matrix"), pybind11::arg("weights"),
                     pybind11::arg("improperRotationIsAllowed") = false,
                     "Initialize a QuaternionFit with a reference, a matrix to "
                     "be fitted and weights");

  quaternion_fit.def(pybind11::init<Eigen::MatrixXd, Eigen::MatrixXd, ElementTypeCollection, bool>(),
                     pybind11::arg("ref_matrix"), pybind11::arg("fit_matrix"), pybind11::arg("elements"),
                     pybind11::arg("improperRotationIsAllowed") = false,
                     "Initialize a QuaternionFit with a reference a to be fitted and use elements as weights");

  quaternion_fit.def("get_rotation_matrix", &QuaternionFit::getRotationMatrix,
                     "Getter for the reverse of the applied "
                     "rotation.");

  quaternion_fit.def("get_trans_vector", &QuaternionFit::getTransVector,
                     "Getter for the reverse of the applied "
                     "translation.");

  quaternion_fit.def("get_fitted_data", &QuaternionFit::getFittedData, "Getter for the fitted data as matrix.");

  quaternion_fit.def("get_rmsd", &QuaternionFit::getRMSD, "Getter for the RMSD not using any weights that might be stored.");

  quaternion_fit.def("get_rot_rmsd", &QuaternionFit::getRotRMSD, "Getter for the RMSD due to differences in rotation only.");

  quaternion_fit.def("get_weighted_rmsd",
                     pybind11::overload_cast<const Eigen::VectorXd&>(&QuaternionFit::getWeightedRMSD, pybind11::const_),
                     "Getter for the RMSD using the given weights.");

  quaternion_fit.def("get_weighted_rmsd", pybind11::overload_cast<>(&QuaternionFit::getWeightedRMSD, pybind11::const_),
                     "Getter for the RMSD using the internal weights given/implied in the constructor.");
}
