/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Math/AtomicSecondDerivativeCollection.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_atomic_second_derivative_collection(pybind11::module& m) {
  pybind11::class_<AtomicSecondDerivativeCollection> atomic_second_derivative_collection(
      m, "AtomicSecondDerivativeCollection");
  atomic_second_derivative_collection.def("get_atomic_hessians", &AtomicSecondDerivativeCollection::getAtomicHessians,
                                          "Get list of atomic hessians");
}
