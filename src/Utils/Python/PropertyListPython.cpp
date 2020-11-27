/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/CalculatorBasics/PropertyList.h>
#include <pybind11/pybind11.h>

using namespace Scine::Utils;

bool property_list_contains(PropertyList& list, Property p) {
  return list.containsSubSet(PropertyList(p));
}

void init_property_list(pybind11::module& m) {
  pybind11::enum_<Property> property(m, "Property");

  property.value("Energy", Property::Energy);
  property.value("Gradients", Property::Gradients);
  property.value("Hessian", Property::Hessian);
  property.value("Dipole", Property::Dipole);
  property.value("DipoleGradient", Property::DipoleGradient);
  property.value("DipoleMatrixAO", Property::DipoleMatrixAO);
  property.value("DipoleMatrixMO", Property::DipoleMatrixMO);
  property.value("DensityMatrix", Property::DensityMatrix);
  property.value("OneElectronMatrix", Property::OneElectronMatrix);
  property.value("TwoElectronMatrix", Property::TwoElectronMatrix);
  property.value("OverlapMatrix", Property::OverlapMatrix);
  property.value("CoefficientMatrix", Property::CoefficientMatrix);
  property.value("BondOrderMatrix", Property::BondOrderMatrix);
  property.value("Thermochemistry", Property::Thermochemistry);
  property.value("ExcitedStates", Property::ExcitedStates);
  property.value("AtomicCharges", Property::AtomicCharges);
  property.value("AtomicGtos", Property::AtomicGtos);
  property.value("AOtoAtomMapping", Property::AOtoAtomMapping);
  property.value("PointChargesGradients", Property::PointChargesGradients);
  property.value("Description", Property::Description);
  property.value("SuccessfulCalculation", Property::SuccessfulCalculation);
  property.value("ProgramName", Property::ProgramName);

  pybind11::class_<PropertyList> property_list(m, "PropertyList");
  property_list.def(pybind11::init<>(), "Empty-initialize");
  property_list.def(pybind11::init<Property>(), "Default initalize from a property");
  property_list.def("add_property", &PropertyList::addProperty, "Add a property to the list");
  property_list.def(
      "contains_subset", &PropertyList::containsSubSet,
      "Check whether this list of properties encompasses all properties set the supplied list of properties");

  property_list.def("__contains__", &property_list_contains);
}
