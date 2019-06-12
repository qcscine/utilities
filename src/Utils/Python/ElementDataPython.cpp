/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/ElementData.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;
using namespace Scine::Utils::Constants;

void init_element_data(pybind11::module& m) {
  pybind11::class_<ElementDataSingleton> element_data_singleton(m, "ElementDataSingleton");

  // Following singleton advice from https://github.com/pybind/pybind11/issues/1418
  element_data_singleton.def_static("instance", &ElementDataSingleton::instance,
                                    pybind11::return_value_policy::reference, "Access the singleton instance");

  // Information access via element type
  element_data_singleton.def("__getitem__",
                             pybind11::overload_cast<const ElementType&>(&ElementDataSingleton::operator[], pybind11::const_),
                             "Access a particular element's data via its element type");

  // Information access via symbol string
  element_data_singleton.def("__getitem__",
                             pybind11::overload_cast<const std::string&>(&ElementDataSingleton::operator[], pybind11::const_),
                             "Access a particular element's data via its symbol string");

  // Nested ElementData class
  pybind11::class_<ElementDataSingleton::ElementData> element_data(element_data_singleton, "ElementData");

  // RO Properties
  element_data.def_property_readonly("symbol", &ElementDataSingleton::ElementData::symbol,
                                     "Two-character string representation of the elements name");
  element_data.def_property_readonly("Z", &ElementDataSingleton::ElementData::Z, "Atomic number of the element");
  element_data.def_property_readonly("mass", &ElementDataSingleton::ElementData::mass, "Relative atomic mass (u)");
  element_data.def_property_readonly("vdw_radius", &ElementDataSingleton::ElementData::vdWRadius,
                                     "van der Waals radius in Angstrom");
  element_data.def_property_readonly("val_electrons", &ElementDataSingleton::ElementData::valElectrons,
                                     "Number of valence electrons");
  element_data.def_property_readonly("s_electrons", &ElementDataSingleton::ElementData::sElectrons,
                                     "Number of valence s electrons");
  element_data.def_property_readonly("p_electrons", &ElementDataSingleton::ElementData::pElectrons,
                                     "Number of valence p electrons");
  element_data.def_property_readonly("d_electrons", &ElementDataSingleton::ElementData::dElectrons,
                                     "Number of valence d electrons");
}
