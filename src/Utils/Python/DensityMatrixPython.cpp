/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/optional.hpp>

using namespace Scine::Utils;

void init_density_matrix(pybind11::module& m) {
  pybind11::class_<DensityMatrix> density_matrix(m, "DensityMatrix");
  density_matrix.def_property_readonly(
      "alpha_matrix", &DensityMatrix::alphaMatrix,
      "Returns the density matrix for electrons with alpha spin from an unrestricted calculation.");
  density_matrix.def_property_readonly(
      "beta_matrix", &DensityMatrix::betaMatrix,
      "Returns the density matrix for electrons with beta spin from an unrestricted calculation.");
  density_matrix.def_property_readonly("restricted_matrix", &DensityMatrix::restrictedMatrix,
                                       "Returns the density matrix for all electrons from a restricted calculation.");
}
