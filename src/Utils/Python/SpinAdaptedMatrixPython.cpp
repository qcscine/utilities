/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/optional.hpp>

using namespace Scine::Utils;

void init_spin_adapted_matrix(pybind11::module& m) {
  pybind11::class_<SpinAdaptedMatrix> spinAdaptedMatrix(m, "SpinAdaptedMatrix");

  spinAdaptedMatrix.def(
      "alpha", [](const SpinAdaptedMatrix& matrix, unsigned i, unsigned j) -> double { return matrix.alpha(i, j); });
  spinAdaptedMatrix.def(
      "beta", [](const SpinAdaptedMatrix& matrix, unsigned i, unsigned j) -> double { return matrix.beta(i, j); });
  spinAdaptedMatrix.def("restricted", [](const SpinAdaptedMatrix& matrix, unsigned i, unsigned j) -> double {
    return matrix.restricted(i, j);
  });
}
