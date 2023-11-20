/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Typenames.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

using namespace Scine::Utils;

void init_typenames(pybind11::module& m) {
  // STL types
  pybind11::bind_vector<ElementTypeCollection>(m, "ElementTypeCollection");
}
