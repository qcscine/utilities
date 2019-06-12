/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/ElementInfo.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_element_info(pybind11::module& m) {
  pybind11::class_<ElementInfo> element_info(m, "ElementInfo");

  element_info.def_static("element_from_symbol", &ElementInfo::elementTypeForSymbol,
                          "Translate a string element type representation to an ElementType");
  element_info.def_static("symbol", &ElementInfo::symbol, "Translate an ElementType into its string representation");
  element_info.def_static("mass", &ElementInfo::mass, "Relative atomic mass of an element");
  element_info.def_static("vdw_radius", &ElementInfo::vdwRadius, "van der Waals radius in Angstrom");
  element_info.def_static("Z", &ElementInfo::Z, "Atomic number");
  element_info.def_static("val_electrons", &ElementInfo::valElectrons, "Number of valence electrons");
  element_info.def_static("s_electrons", &ElementInfo::sElectrons, "Number of s valence electrons");
  element_info.def_static("p_electrons", &ElementInfo::pElectrons, "Number of p valence electrons");
  element_info.def_static("d_electrons", &ElementInfo::dElectrons, "Number of d valence electrons");
}
