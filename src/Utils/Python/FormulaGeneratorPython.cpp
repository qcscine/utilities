/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/FormulaGenerator.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_formula_generator(pybind11::module& m) {
  m.def("generate_chemical_formula", &generateChemicalFormula, pybind11::arg("elements"), pybind11::arg("number_prefix") = "",
        pybind11::arg("number_postfix") = "", "Returns a string of a compound's elemental composition");
}
