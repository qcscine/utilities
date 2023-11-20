/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/DataStructures/PartialHessian.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_partial_hessian(pybind11::module& m) {
  pybind11::class_<PartialHessian> partial_hessian(
      m, "PartialHessian",
      "Class defining a partial Hessian for use in embedding calculations. It is simply a container for the matrix and "
      "the indices of the atoms in the supersystem.");

  partial_hessian.def_property_readonly("matrix", &PartialHessian::getMatrix, "Returns the Hessian matrix.");
  partial_hessian.def_property_readonly(
      "indices", &PartialHessian::getIndices,
      "Returns the indices of the atoms part of the partial system for within the super system.");
}
