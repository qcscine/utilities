/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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
  pybind11::class_<SpinAdaptedMatrix> spin_adapted_matrix(m, "SpinAdaptedMatrix");

  spin_adapted_matrix.def_property_readonly(
      "alpha_matrix", pybind11::overload_cast<>(&SpinAdaptedMatrix::alphaMatrix, pybind11::const_),
      "Returns the two-electron matrix of all electrons with alpha spin from an unrestricted calculation.");
  spin_adapted_matrix.def_property_readonly(
      "beta_matrix", pybind11::overload_cast<>(&SpinAdaptedMatrix::betaMatrix, pybind11::const_),
      "Returns the two-electron matrix of all electrons with beta spin from an unrestricted calculation.");
  spin_adapted_matrix.def_property_readonly(
      "restricted_matrix", pybind11::overload_cast<>(&SpinAdaptedMatrix::restrictedMatrix, pybind11::const_),
      "Returns the two-electron matrix for all electrons from a restricted calculation.");
}
